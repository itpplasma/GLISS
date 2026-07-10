module mercier_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_reconstruction, only: periodic_sixth_order_derivatives, &
        reconstruct_harmonic_grid, reconstruction_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, harmonic_pair_t, &
        radial_grid_half
    use nonuniform_derivative, only: first_derivative_nonuniform
    implicit none
    private

    integer, parameter, public :: mercier_ok = 0
    integer, parameter, public :: mercier_invalid_input = 1
    integer, parameter, public :: mercier_reconstruction_error = 2

    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)
    real(dp), parameter :: mu0 = 2.0_dp * two_pi * 1.0e-7_dp

    type, public :: mercier_result_t
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: d_shear(:)
        real(dp), allocatable :: d_current(:)
        real(dp), allocatable :: d_well(:)
        real(dp), allocatable :: d_geodesic(:)
        real(dp), allocatable :: d_mercier(:)
        real(dp), allocatable :: iota_deviation(:)
        real(dp), allocatable :: boozer_deviation(:)
        real(dp), allocatable :: force_balance_residual(:)
    end type mercier_result_t

    type :: surface_data_t
        real(dp), allocatable :: jacobian(:, :), g_tt(:, :), g_tz(:, :)
        real(dp), allocatable :: g_zz(:, :), b_theta(:, :), b_zeta(:, :)
        real(dp), allocatable :: mod_b(:, :), beta_s(:, :)
        real(dp), allocatable :: area_element(:, :)
    end type surface_data_t

    public :: compute_mercier

contains

    subroutine compute_mercier(equilibrium, n_theta, n_zeta, result, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: n_theta, n_zeta
        type(mercier_result_t), intent(out) :: result
        integer, intent(out) :: info
        type(surface_data_t) :: surface
        type(harmonic_pair_t) :: xhat_s, yhat_s, zhat_s
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp), allocatable :: covariant_theta(:), covariant_zeta(:)
        real(dp), allocatable :: flux_slope(:), volume_slope(:)
        real(dp), allocatable :: covariant_theta_slope(:)
        real(dp), allocatable :: covariant_zeta_slope(:)
        real(dp), allocatable :: pressure_slope(:), iota_slope(:)
        real(dp), allocatable :: flux_curvature(:), volume_curvature(:)
        integer :: ns, i

        info = mercier_invalid_input
        if (n_theta < 8 .or. n_zeta < 8) return
        if (equilibrium%radial_grid /= radial_grid_half) return
        ns = size(equilibrium%s)
        if (ns < 5) return

        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        call differentiate_pair(equilibrium%s, equilibrium%xhat, xhat_s)
        call differentiate_pair(equilibrium%s, equilibrium%yhat, yhat_s)
        call differentiate_pair(equilibrium%s, equilibrium%zhat, zhat_s)

        call allocate_result(ns, result)
        allocate (covariant_theta(ns), covariant_zeta(ns), flux_slope(ns))
        allocate (volume_slope(ns), covariant_theta_slope(ns))
        allocate (covariant_zeta_slope(ns), pressure_slope(ns))
        allocate (iota_slope(ns), flux_curvature(ns), volume_curvature(ns))

        do i = 1, ns
            call load_surface(equilibrium, xhat_s, yhat_s, zhat_s, i, &
                theta, zeta, surface, info)
            if (info /= mercier_ok) return
            covariant_theta(i) = grid_mean(surface%g_tt * surface%b_theta &
                + surface%g_tz * surface%b_zeta)
            covariant_zeta(i) = grid_mean(surface%g_tz * surface%b_theta &
                + surface%g_zz * surface%b_zeta)
            flux_slope(i) = grid_mean(surface%jacobian * surface%b_zeta)
            volume_slope(i) = grid_mean(abs(surface%jacobian)) &
                * real(equilibrium%field_periods, dp)
            result%iota_deviation(i) = abs(real(equilibrium%field_periods, &
                dp) * grid_mean(surface%b_theta / surface%b_zeta) &
                - equilibrium%rotational_transform(i))
            result%boozer_deviation(i) = boozer_deviation(surface, &
                covariant_theta(i), covariant_zeta(i))
        end do

        call first_derivative_nonuniform(equilibrium%s, &
            equilibrium%pressure, pressure_slope)
        call first_derivative_nonuniform(equilibrium%s, &
            equilibrium%rotational_transform, iota_slope)
        call first_derivative_nonuniform(equilibrium%s, covariant_theta, &
            covariant_theta_slope)
        call first_derivative_nonuniform(equilibrium%s, covariant_zeta, &
            covariant_zeta_slope)
        call first_derivative_nonuniform(equilibrium%s, flux_slope, &
            flux_curvature)
        call first_derivative_nonuniform(equilibrium%s, volume_slope, &
            volume_curvature)

        do i = 1, ns
            call load_surface(equilibrium, xhat_s, yhat_s, zhat_s, i, &
                theta, zeta, surface, info)
            if (info /= mercier_ok) return
            call assemble_surface_terms(equilibrium, surface, i, &
                covariant_theta(i), covariant_zeta(i), &
                covariant_theta_slope(i), covariant_zeta_slope(i), &
                flux_slope(i), flux_curvature(i), volume_slope(i), &
                volume_curvature(i), pressure_slope(i), iota_slope(i), &
                result)
        end do
        result%s = equilibrium%s
        info = mercier_ok
    end subroutine compute_mercier

    subroutine assemble_surface_terms(equilibrium, surface, i, &
            covariant_theta, covariant_zeta, covariant_theta_slope, &
            covariant_zeta_slope, flux_slope, flux_curvature, &
            volume_slope, volume_curvature, pressure_slope, iota_slope, &
            result)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(surface_data_t), intent(in) :: surface
        integer, intent(in) :: i
        real(dp), intent(in) :: covariant_theta, covariant_zeta
        real(dp), intent(in) :: covariant_theta_slope, covariant_zeta_slope
        real(dp), intent(in) :: flux_slope, flux_curvature
        real(dp), intent(in) :: volume_slope, volume_curvature
        real(dp), intent(in) :: pressure_slope, iota_slope
        type(mercier_result_t), intent(inout) :: result
        real(dp), allocatable :: beta_theta(:, :), beta_zeta(:, :)
        real(dp), allocatable :: mu0_j_dot_b(:, :), grad_psi(:, :)
        real(dp), allocatable :: b_squared(:, :)
        real(dp) :: field_periods, psi_slope, psi_curvature, iota_psi
        real(dp) :: dp_dpsi, d2v_dpsi2, current_slope_ratio
        real(dp) :: integral_xi, integral_inverse, integral_bsq
        real(dp) :: integral_jb, integral_jb_squared, n_grid
        real(dp) :: pressure_term, toroidal_term, poloidal_term
        integer :: fd_info

        field_periods = real(equilibrium%field_periods, dp)
        n_grid = real(size(surface%jacobian), dp)
        call periodic_sixth_order_derivatives(surface%beta_s, &
            1.0_dp / real(size(surface%beta_s, 1), dp), &
            1.0_dp / real(size(surface%beta_s, 2), dp), beta_theta, &
            beta_zeta, fd_info)
        mu0_j_dot_b = ((beta_zeta - covariant_zeta_slope) &
            * covariant_theta &
            + (covariant_theta_slope - beta_theta) * covariant_zeta) &
            / surface%jacobian
        b_squared = surface%mod_b**2
        psi_slope = flux_slope / two_pi
        psi_curvature = flux_curvature / two_pi
        iota_psi = iota_slope / psi_slope
        dp_dpsi = mu0 * pressure_slope / psi_slope
        grad_psi = abs(psi_slope) * surface%area_element &
            / abs(surface%jacobian)

        result%d_shear(i) = iota_psi**2 / (16.0_dp * acos(-1.0_dp)**2)
        current_slope_ratio = covariant_theta_slope / (two_pi * psi_slope)
        integral_xi = field_periods * sum(surface%area_element &
            * (mu0_j_dot_b - current_slope_ratio * b_squared) &
            / grad_psi**3) / n_grid
        result%d_current(i) = -sign(1.0_dp, covariant_zeta) &
            / two_pi**4 * iota_psi * integral_xi
        d2v_dpsi2 = (volume_curvature * psi_slope &
            - volume_slope * psi_curvature) / psi_slope**3
        integral_inverse = field_periods * sum(surface%area_element &
            / (b_squared * grad_psi)) / n_grid
        integral_bsq = field_periods * sum(surface%area_element &
            * b_squared / grad_psi**3) / n_grid
        result%d_well(i) = dp_dpsi * (sign(1.0_dp, psi_slope) * d2v_dpsi2 &
            - dp_dpsi * integral_inverse) * integral_bsq / two_pi**6
        integral_jb = field_periods * sum(surface%area_element &
            * mu0_j_dot_b / grad_psi**3) / n_grid
        integral_jb_squared = field_periods * sum(surface%area_element &
            * mu0_j_dot_b**2 / (b_squared * grad_psi**3)) / n_grid
        result%d_geodesic(i) = (integral_jb**2 &
            - integral_bsq * integral_jb_squared) / two_pi**6
        result%d_mercier(i) = result%d_shear(i) + result%d_current(i) &
            + result%d_well(i) + result%d_geodesic(i)
        pressure_term = mu0 * pressure_slope * grid_mean(surface%jacobian)
        toroidal_term = flux_slope * covariant_zeta_slope
        poloidal_term = grid_mean(surface%jacobian * surface%b_theta) &
            * covariant_theta_slope
        result%force_balance_residual(i) = &
            abs(pressure_term + toroidal_term + poloidal_term) &
            / max(abs(pressure_term) + abs(toroidal_term) &
            + abs(poloidal_term), &
            1.0e-14_dp * abs(flux_slope * covariant_zeta))
    end subroutine assemble_surface_terms

    subroutine load_surface(equilibrium, xhat_s, yhat_s, zhat_s, i, theta, &
            zeta, surface, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(harmonic_pair_t), intent(in) :: xhat_s, yhat_s, zhat_s
        integer, intent(in) :: i
        real(dp), intent(in) :: theta(:), zeta(:)
        type(surface_data_t), intent(out) :: surface
        integer, intent(out) :: info
        real(dp), allocatable :: es_x(:, :), es_y(:, :), es_z(:, :)
        real(dp), allocatable :: eu_x(:, :), eu_y(:, :), eu_z(:, :)
        real(dp), allocatable :: ev_x(:, :), ev_y(:, :), ev_z(:, :)
        real(dp), allocatable :: discard_a(:, :), discard_b(:, :)
        real(dp), allocatable :: g_su(:, :), g_sv(:, :)

        call surface_values(equilibrium%jacobian, i, equilibrium, theta, &
            zeta, surface%jacobian, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%g_tt, i, equilibrium, theta, zeta, &
            surface%g_tt, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%g_tz, i, equilibrium, theta, zeta, &
            surface%g_tz, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%g_zz, i, equilibrium, theta, zeta, &
            surface%g_zz, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%b_contravariant_theta, i, &
            equilibrium, theta, zeta, surface%b_theta, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%b_contravariant_zeta, i, &
            equilibrium, theta, zeta, surface%b_zeta, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%mod_b, i, equilibrium, theta, zeta, &
            surface%mod_b, info)
        if (info /= mercier_ok) return

        call surface_with_derivatives(xhat_s, i, equilibrium, theta, zeta, &
            es_x, discard_a, discard_b, info)
        if (info /= mercier_ok) return
        call surface_with_derivatives(yhat_s, i, equilibrium, theta, zeta, &
            es_y, discard_a, discard_b, info)
        if (info /= mercier_ok) return
        call surface_with_derivatives(zhat_s, i, equilibrium, theta, zeta, &
            es_z, discard_a, discard_b, info)
        if (info /= mercier_ok) return
        call surface_with_derivatives(equilibrium%xhat, i, equilibrium, &
            theta, zeta, discard_a, eu_x, ev_x, info)
        if (info /= mercier_ok) return
        call surface_with_derivatives(equilibrium%yhat, i, equilibrium, &
            theta, zeta, discard_a, eu_y, ev_y, info)
        if (info /= mercier_ok) return
        call surface_with_derivatives(equilibrium%zhat, i, equilibrium, &
            theta, zeta, discard_a, eu_z, ev_z, info)
        if (info /= mercier_ok) return

        g_su = es_x * eu_x + es_y * eu_y + es_z * eu_z
        g_sv = es_x * ev_x + es_y * ev_y + es_z * ev_z
        surface%beta_s = g_su * surface%b_theta + g_sv * surface%b_zeta
        surface%area_element = sqrt(max(surface%g_tt * surface%g_zz &
            - surface%g_tz**2, 0.0_dp))
        info = mercier_ok
    end subroutine load_surface

    subroutine surface_values(pair, i, equilibrium, theta, zeta, values, &
            info)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: i
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), allocatable, intent(out) :: values(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: discard_theta(:, :), discard_zeta(:, :)
        integer :: rec_info

        call reconstruct_harmonic_grid(pair, i, equilibrium%poloidal_modes, &
            equilibrium%toroidal_modes, theta, zeta, values, discard_theta, &
            discard_zeta, rec_info)
        info = merge(mercier_ok, mercier_reconstruction_error, &
            rec_info == reconstruction_ok)
    end subroutine surface_values

    subroutine surface_with_derivatives(pair, i, equilibrium, theta, zeta, &
            values, derivative_theta, derivative_zeta, info)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: i
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), allocatable, intent(out) :: values(:, :)
        real(dp), allocatable, intent(out) :: derivative_theta(:, :)
        real(dp), allocatable, intent(out) :: derivative_zeta(:, :)
        integer, intent(out) :: info
        integer :: rec_info

        call reconstruct_harmonic_grid(pair, i, equilibrium%poloidal_modes, &
            equilibrium%toroidal_modes, theta, zeta, values, &
            derivative_theta, derivative_zeta, rec_info)
        info = merge(mercier_ok, mercier_reconstruction_error, &
            rec_info == reconstruction_ok)
    end subroutine surface_with_derivatives

    pure subroutine differentiate_pair(s, pair, slope_pair)
        real(dp), intent(in) :: s(:)
        type(harmonic_pair_t), intent(in) :: pair
        type(harmonic_pair_t), intent(out) :: slope_pair
        integer :: m, n

        allocate (slope_pair%cosine, mold=pair%cosine)
        allocate (slope_pair%sine, mold=pair%sine)
        do n = 1, size(pair%cosine, 3)
            do m = 1, size(pair%cosine, 2)
                call first_derivative_nonuniform(s, pair%cosine(:, m, n), &
                    slope_pair%cosine(:, m, n))
                call first_derivative_nonuniform(s, pair%sine(:, m, n), &
                    slope_pair%sine(:, m, n))
            end do
        end do
    end subroutine differentiate_pair

    pure function boozer_deviation(surface, covariant_theta, &
            covariant_zeta) result(deviation)
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: covariant_theta, covariant_zeta
        real(dp) :: deviation
        real(dp) :: scale

        scale = max(abs(covariant_theta), abs(covariant_zeta))
        deviation = max(maxval(abs(surface%g_tt * surface%b_theta &
            + surface%g_tz * surface%b_zeta - covariant_theta)), &
            maxval(abs(surface%g_tz * surface%b_theta &
            + surface%g_zz * surface%b_zeta - covariant_zeta))) / scale
    end function boozer_deviation

    pure function grid_mean(values) result(mean)
        real(dp), intent(in) :: values(:, :)
        real(dp) :: mean

        mean = sum(values) / real(size(values), dp)
    end function grid_mean

    pure subroutine build_angular_grids(n_theta, n_zeta, theta, zeta)
        integer, intent(in) :: n_theta, n_zeta
        real(dp), allocatable, intent(out) :: theta(:), zeta(:)
        integer :: j

        allocate (theta(n_theta), zeta(n_zeta))
        do j = 1, n_theta
            theta(j) = real(j - 1, dp) / real(n_theta, dp)
        end do
        do j = 1, n_zeta
            zeta(j) = real(j - 1, dp) / real(n_zeta, dp)
        end do
    end subroutine build_angular_grids

    pure subroutine allocate_result(ns, result)
        integer, intent(in) :: ns
        type(mercier_result_t), intent(out) :: result

        allocate (result%s(ns), result%d_shear(ns), result%d_current(ns))
        allocate (result%d_well(ns), result%d_geodesic(ns))
        allocate (result%d_mercier(ns), result%iota_deviation(ns))
        allocate (result%boozer_deviation(ns))
        allocate (result%force_balance_residual(ns))
    end subroutine allocate_result

end module mercier_diagnostic
