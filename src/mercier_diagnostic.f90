module mercier_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_reconstruction, only: project_harmonic_grid, &
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

    type, public :: mercier_d_terms_t
        real(dp) :: shear, current, well, geodesic, mercier
    end type mercier_d_terms_t

    type, public :: mercier_gradient_t
        real(dp) :: d_iota_slope, d_pressure_slope
    end type mercier_gradient_t

    type, public :: mercier_result_t
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: d_shear(:)
        real(dp), allocatable :: d_current(:)
        real(dp), allocatable :: d_well(:)
        real(dp), allocatable :: d_geodesic(:)
        real(dp), allocatable :: d_mercier(:)
        real(dp), allocatable :: d_mercier_d_iota_slope(:)
        real(dp), allocatable :: d_mercier_d_pressure_slope(:)
        real(dp), allocatable :: d_mercier_d_pressure_slope_full(:)
        real(dp), allocatable :: d_mercier_d_toroidal_flux_slope(:)
        real(dp), allocatable :: d_mercier_d_poloidal_flux_slope(:)
        real(dp), allocatable :: iota_deviation(:)
        real(dp), allocatable :: boozer_deviation(:)
        real(dp), allocatable :: force_balance_residual(:)
        real(dp), allocatable :: jacobian_identity_deviation(:)
        real(dp), allocatable :: beta_chart_deviation(:)
    end type mercier_result_t

    type, public :: surface_data_t
        real(dp), allocatable :: jacobian(:, :), g_tt(:, :), g_tz(:, :)
        real(dp), allocatable :: g_zz(:, :), b_theta(:, :), b_zeta(:, :)
        real(dp), allocatable :: g_st(:, :), g_sz(:, :)
        real(dp), allocatable :: mod_b(:, :)
        real(dp), allocatable :: area_element(:, :)
    end type surface_data_t

    type :: surface_profiles_t
        real(dp) :: flux_slope, poloidal_slope
        real(dp) :: flux_curvature, poloidal_curvature
        real(dp) :: covariant_theta, covariant_zeta
        real(dp) :: covariant_theta_slope, covariant_zeta_slope
        real(dp) :: pressure_slope
    end type surface_profiles_t

    public :: compute_mercier
    public :: build_angular_grids
    public :: build_kernel_geometry
    public :: differentiate_pair
    public :: mercier_d_terms
    public :: mercier_d_terms_gradient
    public :: mercier_surface_terms

contains

    subroutine build_kernel_geometry(equilibrium, n_theta, n_zeta, &
            fields, drive, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: n_theta, n_zeta
        real(dp), allocatable, intent(out) :: fields(:, :, :, :)
        real(dp), allocatable, intent(out) :: drive(:, :, :)
        integer, intent(out) :: info
        type(surface_data_t) :: surface
        type(harmonic_pair_t) :: jacobian_slope
        type(surface_profiles_t) :: profiles
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp), allocatable :: covariant_theta(:), covariant_zeta(:)
        real(dp), allocatable :: flux_slope(:), poloidal_slope(:)
        real(dp), allocatable :: covariant_theta_slope(:)
        real(dp), allocatable :: covariant_zeta_slope(:)
        real(dp), allocatable :: pressure_slope(:)
        real(dp), allocatable :: flux_curvature(:), poloidal_curvature(:)
        integer :: ns, i

        info = mercier_invalid_input
        if (n_theta < 8 .or. n_zeta < 8) return
        if (equilibrium%radial_grid /= radial_grid_half) return
        ns = size(equilibrium%s)
        if (ns < 5) return

        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        call differentiate_pair(equilibrium%s, equilibrium%jacobian, &
            jacobian_slope)
        allocate (covariant_theta(ns), covariant_zeta(ns))
        allocate (flux_slope(ns), poloidal_slope(ns))
        allocate (covariant_theta_slope(ns), covariant_zeta_slope(ns))
        allocate (pressure_slope(ns), flux_curvature(ns))
        allocate (poloidal_curvature(ns))
        allocate (fields(n_theta, n_zeta, 13, ns))
        allocate (drive(n_theta, n_zeta, ns))

        do i = 1, ns
            call load_surface(equilibrium, i, theta, zeta, surface, info)
            if (info /= mercier_ok) return
            covariant_theta(i) = grid_mean(surface%g_tt * surface%b_theta &
                + surface%g_tz * surface%b_zeta)
            covariant_zeta(i) = grid_mean(surface%g_tz * surface%b_theta &
                + surface%g_zz * surface%b_zeta)
            flux_slope(i) = grid_mean(surface%jacobian * surface%b_zeta)
            poloidal_slope(i) = grid_mean(surface%jacobian &
                * surface%b_theta)
        end do
        call first_derivative_nonuniform(equilibrium%s, &
            equilibrium%pressure, pressure_slope)
        call first_derivative_nonuniform(equilibrium%s, covariant_theta, &
            covariant_theta_slope)
        call first_derivative_nonuniform(equilibrium%s, covariant_zeta, &
            covariant_zeta_slope)
        call first_derivative_nonuniform(equilibrium%s, flux_slope, &
            flux_curvature)
        call first_derivative_nonuniform(equilibrium%s, poloidal_slope, &
            poloidal_curvature)

        do i = 1, ns
            call load_surface(equilibrium, i, theta, zeta, surface, info)
            if (info /= mercier_ok) return
            profiles = surface_profiles_t(flux_slope(i), &
                poloidal_slope(i), flux_curvature(i), &
                poloidal_curvature(i), covariant_theta(i), &
                covariant_zeta(i), covariant_theta_slope(i), &
                covariant_zeta_slope(i), pressure_slope(i))
            call fill_surface_fields(equilibrium, surface, profiles, &
                jacobian_slope, i, theta, zeta, fields(:, :, :, i), &
                drive(:, :, i), info)
            if (info /= mercier_ok) return
        end do
        info = mercier_ok
    end subroutine build_kernel_geometry

    subroutine fill_surface_fields(equilibrium, surface, profiles, &
            jacobian_slope, i, theta, zeta, fields, drive, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(surface_data_t), intent(in) :: surface
        type(surface_profiles_t), intent(in) :: profiles
        type(harmonic_pair_t), intent(in) :: jacobian_slope
        integer, intent(in) :: i
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(out) :: fields(:, :, :)
        real(dp), intent(out) :: drive(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: beta_values(:, :), beta_theta(:, :)
        real(dp), allocatable :: beta_zeta(:, :), jac_slope_grid(:, :)
        real(dp), allocatable :: discard_a(:, :), discard_b(:, :)
        integer :: rec_info

        call solve_beta_derivatives(equilibrium, surface, theta, zeta, &
            profiles%covariant_theta_slope, &
            profiles%covariant_zeta_slope, profiles%pressure_slope, &
            profiles%poloidal_slope, profiles%flux_slope, beta_values, &
            beta_theta, beta_zeta)
        call reconstruct_harmonic_grid(jacobian_slope, i, &
            equilibrium%poloidal_modes, equilibrium%toroidal_modes, &
            theta, zeta, jac_slope_grid, discard_a, discard_b, rec_info)
        if (rec_info /= reconstruction_ok) then
            info = mercier_reconstruction_error
            return
        end if
        fields(:, :, 1) = profiles%flux_slope
        fields(:, :, 2) = profiles%poloidal_slope
        fields(:, :, 3) = profiles%flux_curvature
        fields(:, :, 4) = profiles%poloidal_curvature
        fields(:, :, 5) = profiles%covariant_zeta
        fields(:, :, 6) = profiles%covariant_theta
        fields(:, :, 7) = surface%jacobian
        fields(:, :, 8) = surface%mod_b
        fields(:, :, 9) = max((surface%g_tt * surface%g_zz &
            - surface%g_tz**2) / surface%jacobian**2, tiny(1.0_dp))
        fields(:, :, 10) = ((beta_zeta &
            - profiles%covariant_zeta_slope) * profiles%covariant_theta &
            + (profiles%covariant_theta_slope - beta_theta) &
            * profiles%covariant_zeta) / surface%jacobian
        fields(:, :, 11) = mu0 * profiles%pressure_slope
        if (equilibrium%has_chart_metric) then
            fields(:, :, 12) = ((surface%g_tz * surface%b_theta &
                + surface%g_zz * surface%b_zeta) * surface%g_st &
                - (surface%g_tt * surface%b_theta &
                + surface%g_tz * surface%b_zeta) * surface%g_sz) &
                / (surface%jacobian * surface%mod_b)
            fields(:, :, 13) = surface%b_theta * surface%g_st &
                + surface%b_zeta * surface%g_sz
        else
            fields(:, :, 12) = 0.0_dp
            fields(:, :, 13) = beta_values
        end if
        ! The current-curvature group enters with a minus sign in the
        ! export chart; pinned against the geometric drive in
        ! derivations/drive_machinery_identity.wl.
        drive = (fields(:, :, 10)**2 &
            + (mu0 * profiles%pressure_slope)**2 * fields(:, :, 9)) &
            / (surface%mod_b**2 * fields(:, :, 9)) &
            - (profiles%flux_curvature * profiles%covariant_zeta_slope &
            + profiles%poloidal_curvature &
            * profiles%covariant_theta_slope &
            - profiles%flux_curvature * beta_zeta &
            - profiles%poloidal_curvature * beta_theta) &
            / surface%jacobian &
            - mu0 * profiles%pressure_slope * jac_slope_grid &
            / surface%jacobian
        info = mercier_ok
        if (equilibrium%has_chart_metric) then
            call add_drive_chart_term(equilibrium, surface, theta, zeta, &
                profiles%covariant_theta_slope, &
                profiles%covariant_zeta_slope, beta_theta, beta_zeta, &
                profiles%poloidal_slope, profiles%flux_slope, drive, &
                info)
        end if
    end subroutine fill_surface_fields

    subroutine add_drive_chart_term(equilibrium, surface, theta, zeta, &
            covariant_theta_slope, covariant_zeta_slope, beta_theta, &
            beta_zeta, poloidal_flux_slope, toroidal_flux_slope, drive, &
            info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(in) :: covariant_theta_slope, covariant_zeta_slope
        real(dp), intent(in) :: beta_theta(:, :), beta_zeta(:, :)
        real(dp), intent(in) :: poloidal_flux_slope, toroidal_flux_slope
        real(dp), intent(inout) :: drive(:, :)
        integer, intent(out) :: info
        type(harmonic_pair_t) :: operand_pair
        real(dp), allocatable :: operand(:, :), grad_s2(:, :)
        real(dp), allocatable :: upper_ts(:, :), upper_zs(:, :)
        real(dp), allocatable :: operand_values(:, :)
        real(dp), allocatable :: operand_theta(:, :), operand_zeta(:, :)
        real(dp) :: operand_cosine(size(equilibrium%poloidal_modes), &
            size(equilibrium%toroidal_modes))
        real(dp) :: operand_sine(size(equilibrium%poloidal_modes), &
            size(equilibrium%toroidal_modes))
        integer :: rec_info

        upper_ts = (surface%g_sz * surface%g_tz &
            - surface%g_st * surface%g_zz) / surface%jacobian**2
        upper_zs = (surface%g_st * surface%g_tz &
            - surface%g_sz * surface%g_tt) / surface%jacobian**2
        grad_s2 = max((surface%g_tt * surface%g_zz - surface%g_tz**2) &
            / surface%jacobian**2, tiny(1.0_dp))
        operand = ((covariant_theta_slope - beta_theta) * upper_ts &
            - (beta_zeta - covariant_zeta_slope) * upper_zs) / grad_s2
        call project_harmonic_grid(operand, equilibrium%poloidal_modes, &
            equilibrium%toroidal_modes, theta, zeta, operand_cosine, &
            operand_sine)
        allocate (operand_pair%cosine(1, size(operand_cosine, 1), &
            size(operand_cosine, 2)))
        allocate (operand_pair%sine(1, size(operand_sine, 1), &
            size(operand_sine, 2)))
        operand_pair%cosine(1, :, :) = operand_cosine
        operand_pair%sine(1, :, :) = operand_sine
        call reconstruct_harmonic_grid(operand_pair, 1, &
            equilibrium%poloidal_modes, equilibrium%toroidal_modes, &
            theta, zeta, operand_values, operand_theta, operand_zeta, &
            rec_info)
        if (rec_info /= reconstruction_ok) then
            info = mercier_reconstruction_error
            return
        end if
        drive = drive + (poloidal_flux_slope * operand_theta &
            + toroidal_flux_slope * operand_zeta) / surface%jacobian
        info = mercier_ok
    end subroutine add_drive_chart_term

    subroutine compute_mercier(equilibrium, n_theta, n_zeta, result, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: n_theta, n_zeta
        type(mercier_result_t), intent(out) :: result
        integer, intent(out) :: info
        type(surface_data_t) :: surface
        type(harmonic_pair_t) :: xhat_s, yhat_s, zhat_s
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp), allocatable :: beta_positions(:, :)
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
            call load_surface(equilibrium, i, theta, zeta, surface, info)
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
            result%jacobian_identity_deviation(i) = maxval(abs( &
                surface%mod_b**2 * surface%jacobian &
                - (flux_slope(i) * covariant_zeta(i) &
                + grid_mean(surface%jacobian * surface%b_theta) &
                * covariant_theta(i)))) &
                / abs(flux_slope(i) * covariant_zeta(i))
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
            call load_surface(equilibrium, i, theta, zeta, surface, info)
            if (info /= mercier_ok) return
            call beta_from_positions(equilibrium, xhat_s, yhat_s, &
                zhat_s, i, theta, zeta, surface, beta_positions, info)
            if (info /= mercier_ok) return
            call assemble_surface_terms(equilibrium, surface, theta, &
                zeta, i, &
                covariant_theta(i), covariant_zeta(i), &
                covariant_theta_slope(i), covariant_zeta_slope(i), &
                flux_slope(i), flux_curvature(i), volume_slope(i), &
                volume_curvature(i), pressure_slope(i), iota_slope(i), &
                beta_positions, result)
        end do
        result%s = equilibrium%s
        info = mercier_ok
    end subroutine compute_mercier

    subroutine assemble_surface_terms(equilibrium, surface, theta, zeta, &
            i, covariant_theta, covariant_zeta, covariant_theta_slope, &
            covariant_zeta_slope, flux_slope, flux_curvature, &
            volume_slope, volume_curvature, pressure_slope, iota_slope, &
            beta_positions, result)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(in) :: beta_positions(:, :)
        integer, intent(in) :: i
        real(dp), intent(in) :: covariant_theta, covariant_zeta
        real(dp), intent(in) :: covariant_theta_slope, covariant_zeta_slope
        real(dp), intent(in) :: flux_slope, flux_curvature
        real(dp), intent(in) :: volume_slope, volume_curvature
        real(dp), intent(in) :: pressure_slope, iota_slope
        type(mercier_result_t), intent(inout) :: result
        real(dp), allocatable :: beta_values(:, :), beta_filtered(:, :)
        real(dp) :: pressure_term, toroidal_term, poloidal_term
        real(dp) :: d_pressure_slope_full
        real(dp) :: d_toroidal_flux_slope, d_poloidal_flux_slope
        type(mercier_d_terms_t) :: d_terms
        type(mercier_gradient_t) :: grad

        call mercier_surface_terms(equilibrium, surface, theta, zeta, &
            covariant_theta, covariant_zeta, covariant_theta_slope, &
            covariant_zeta_slope, flux_slope, flux_curvature, &
            volume_slope, volume_curvature, pressure_slope, iota_slope, &
            d_terms, grad, d_pressure_slope_full, beta_values, &
            d_toroidal_flux_slope, d_poloidal_flux_slope)
        result%d_shear(i) = d_terms%shear
        result%d_current(i) = d_terms%current
        result%d_well(i) = d_terms%well
        result%d_geodesic(i) = d_terms%geodesic
        result%d_mercier(i) = d_terms%mercier
        result%d_mercier_d_iota_slope(i) = grad%d_iota_slope
        result%d_mercier_d_pressure_slope(i) = grad%d_pressure_slope
        result%d_mercier_d_pressure_slope_full(i) = d_pressure_slope_full
        result%d_mercier_d_toroidal_flux_slope(i) = d_toroidal_flux_slope
        result%d_mercier_d_poloidal_flux_slope(i) = d_poloidal_flux_slope
        call filter_gauge_kernel(equilibrium, theta, zeta, &
            grid_mean(surface%jacobian * surface%b_theta), flux_slope, &
            beta_positions, beta_filtered)
        result%beta_chart_deviation(i) = &
            maxval(abs(beta_filtered - beta_values)) &
            / max(maxval(abs(beta_values)), &
            1.0e-9_dp * abs(covariant_zeta))
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

    subroutine mercier_surface_terms(equilibrium, surface, theta, zeta, &
            covariant_theta, covariant_zeta, covariant_theta_slope, &
            covariant_zeta_slope, flux_slope, flux_curvature, &
            volume_slope, volume_curvature, pressure_slope, iota_slope, &
            d_terms, grad, d_pressure_slope_full, beta_values, &
            d_toroidal_flux_slope, d_poloidal_flux_slope, &
            poloidal_flux_slope_override)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(in) :: covariant_theta, covariant_zeta
        real(dp), intent(in) :: covariant_theta_slope, covariant_zeta_slope
        real(dp), intent(in) :: flux_slope, flux_curvature
        real(dp), intent(in) :: volume_slope, volume_curvature
        real(dp), intent(in) :: pressure_slope, iota_slope
        type(mercier_d_terms_t), intent(out) :: d_terms
        type(mercier_gradient_t), intent(out) :: grad
        real(dp), intent(out) :: d_pressure_slope_full
        real(dp), allocatable, intent(out) :: beta_values(:, :)
        real(dp), intent(out) :: d_toroidal_flux_slope, d_poloidal_flux_slope
        real(dp), intent(in), optional :: poloidal_flux_slope_override
        real(dp), allocatable :: beta_theta(:, :), beta_zeta(:, :)
        real(dp), allocatable :: mu0_j_dot_b(:, :), grad_psi(:, :)
        real(dp), allocatable :: b_squared(:, :)
        type(harmonic_pair_t) :: beta_harmonics
        real(dp) :: field_periods, psi_slope, psi_curvature
        real(dp) :: poloidal_flux_slope, d2v_dpsi2, current_slope_ratio
        real(dp) :: integral_xi, integral_inverse, integral_bsq
        real(dp) :: integral_jb, integral_jb_squared, n_grid
        real(dp) :: d_integral_mu0jb_toroidal, d_integral_jbsq_toroidal
        real(dp) :: d_integral_mu0jb_poloidal, d_integral_jbsq_poloidal

        field_periods = real(equilibrium%field_periods, dp)
        n_grid = real(size(surface%jacobian), dp)
        poloidal_flux_slope = grid_mean(surface%jacobian * surface%b_theta)
        if (present(poloidal_flux_slope_override)) &
            poloidal_flux_slope = poloidal_flux_slope_override
        call solve_beta_derivatives(equilibrium, surface, theta, zeta, &
            covariant_theta_slope, covariant_zeta_slope, pressure_slope, &
            poloidal_flux_slope, flux_slope, &
            beta_values, beta_theta, beta_zeta, beta_harmonics)
        mu0_j_dot_b = ((beta_zeta - covariant_zeta_slope) &
            * covariant_theta &
            + (covariant_theta_slope - beta_theta) * covariant_zeta) &
            / surface%jacobian
        b_squared = surface%mod_b**2
        psi_slope = flux_slope / two_pi
        psi_curvature = flux_curvature / two_pi
        grad_psi = abs(psi_slope) * surface%area_element &
            / abs(surface%jacobian)

        current_slope_ratio = covariant_theta_slope / (two_pi * psi_slope)
        integral_xi = field_periods * sum(surface%area_element &
            * (mu0_j_dot_b - current_slope_ratio * b_squared) &
            / grad_psi**3) / n_grid
        d2v_dpsi2 = (volume_curvature * psi_slope &
            - volume_slope * psi_curvature) / psi_slope**3
        integral_inverse = field_periods * sum(surface%area_element &
            / (b_squared * grad_psi)) / n_grid
        integral_bsq = field_periods * sum(surface%area_element &
            * b_squared / grad_psi**3) / n_grid
        integral_jb = field_periods * sum(surface%area_element &
            * mu0_j_dot_b / grad_psi**3) / n_grid
        integral_jb_squared = field_periods * sum(surface%area_element &
            * mu0_j_dot_b**2 / (b_squared * grad_psi**3)) / n_grid
        d_terms = mercier_d_terms(iota_slope, pressure_slope, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            integral_bsq, integral_jb, integral_jb_squared)
        grad = mercier_d_terms_gradient(iota_slope, pressure_slope, &
            psi_slope, covariant_zeta, integral_xi, d2v_dpsi2, &
            integral_inverse, integral_bsq)
        d_pressure_slope_full = grad%d_pressure_slope &
            + pressure_implicit_gradient(equilibrium, surface, theta, &
            zeta, poloidal_flux_slope, flux_slope, covariant_theta, &
            covariant_zeta, mu0_j_dot_b, psi_slope, field_periods, &
            n_grid, iota_slope, integral_jb, integral_bsq)
        call flux_beta_integral_derivatives(equilibrium, surface, theta, &
            zeta, beta_harmonics, poloidal_flux_slope, flux_slope, &
            covariant_theta, covariant_zeta, mu0_j_dot_b, grad_psi, &
            b_squared, field_periods, n_grid, d_integral_mu0jb_toroidal, &
            d_integral_jbsq_toroidal, d_integral_mu0jb_poloidal, &
            d_integral_jbsq_poloidal)
        call mercier_flux_slope_gradients(iota_slope, pressure_slope, &
            psi_slope, covariant_zeta, flux_slope, current_slope_ratio, &
            volume_curvature, d2v_dpsi2, integral_xi, integral_inverse, &
            integral_bsq, integral_jb, integral_jb_squared, &
            d_integral_mu0jb_toroidal, d_integral_jbsq_toroidal, &
            d_integral_mu0jb_poloidal, d_integral_jbsq_poloidal, &
            d_toroidal_flux_slope, d_poloidal_flux_slope)
    end subroutine mercier_surface_terms

    function pressure_implicit_gradient(equilibrium, surface, theta, zeta, &
            poloidal_flux_slope, toroidal_flux_slope, covariant_theta, &
            covariant_zeta, mu0_j_dot_b, psi_slope, field_periods, &
            n_grid, iota_slope, integral_jb, integral_bsq) &
            result(contribution)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(in) :: poloidal_flux_slope, toroidal_flux_slope
        real(dp), intent(in) :: covariant_theta, covariant_zeta
        real(dp), intent(in) :: mu0_j_dot_b(:, :), psi_slope
        real(dp), intent(in) :: field_periods, n_grid, iota_slope
        real(dp), intent(in) :: integral_jb, integral_bsq
        real(dp) :: contribution
        real(dp), allocatable :: discard_values(:, :), d_beta_theta(:, :)
        real(dp), allocatable :: d_beta_zeta(:, :), d_mu0_j_dot_b(:, :)
        real(dp), allocatable :: grad_psi(:, :), b_squared(:, :)
        real(dp) :: d_integral_mu0jb, d_integral_jb_squared, iota_psi

        call solve_beta_derivatives(equilibrium, surface, theta, zeta, &
            0.0_dp, 0.0_dp, 1.0_dp, poloidal_flux_slope, &
            toroidal_flux_slope, discard_values, d_beta_theta, d_beta_zeta)
        d_mu0_j_dot_b = (d_beta_zeta * covariant_theta &
            - d_beta_theta * covariant_zeta) / surface%jacobian
        grad_psi = abs(psi_slope) * surface%area_element &
            / abs(surface%jacobian)
        b_squared = surface%mod_b**2
        d_integral_mu0jb = field_periods * sum(surface%area_element &
            * d_mu0_j_dot_b / grad_psi**3) / n_grid
        d_integral_jb_squared = field_periods * sum(surface%area_element &
            * 2.0_dp * mu0_j_dot_b * d_mu0_j_dot_b &
            / (b_squared * grad_psi**3)) / n_grid
        iota_psi = iota_slope / psi_slope
        contribution = (-sign(1.0_dp, covariant_zeta) / two_pi**4 &
            * iota_psi + 2.0_dp * integral_jb / two_pi**6) &
            * d_integral_mu0jb &
            - integral_bsq / two_pi**6 * d_integral_jb_squared
    end function pressure_implicit_gradient

    subroutine flux_beta_integral_derivatives(equilibrium, surface, theta, &
            zeta, beta_harmonics, poloidal_flux_slope, toroidal_flux_slope, &
            covariant_theta, covariant_zeta, mu0_j_dot_b, grad_psi, &
            b_squared, field_periods, n_grid, d_integral_mu0jb_toroidal, &
            d_integral_jbsq_toroidal, d_integral_mu0jb_poloidal, &
            d_integral_jbsq_poloidal)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: theta(:), zeta(:)
        type(harmonic_pair_t), intent(in) :: beta_harmonics
        real(dp), intent(in) :: poloidal_flux_slope, toroidal_flux_slope
        real(dp), intent(in) :: covariant_theta, covariant_zeta
        real(dp), intent(in) :: mu0_j_dot_b(:, :), grad_psi(:, :)
        real(dp), intent(in) :: b_squared(:, :), field_periods, n_grid
        real(dp), intent(out) :: d_integral_mu0jb_toroidal
        real(dp), intent(out) :: d_integral_jbsq_toroidal
        real(dp), intent(out) :: d_integral_mu0jb_poloidal
        real(dp), intent(out) :: d_integral_jbsq_poloidal
        real(dp), allocatable :: d_beta_theta(:, :), d_beta_zeta(:, :)
        real(dp), allocatable :: d_mu0_toroidal(:, :), d_mu0_poloidal(:, :)

        call beta_flux_slope_derivative(equilibrium, beta_harmonics, &
            poloidal_flux_slope, toroidal_flux_slope, 0.0_dp, 1.0_dp, &
            theta, zeta, d_beta_theta, d_beta_zeta)
        d_mu0_toroidal = (d_beta_zeta * covariant_theta &
            - d_beta_theta * covariant_zeta) / surface%jacobian
        call beta_flux_slope_derivative(equilibrium, beta_harmonics, &
            poloidal_flux_slope, toroidal_flux_slope, 1.0_dp, 0.0_dp, &
            theta, zeta, d_beta_theta, d_beta_zeta)
        d_mu0_poloidal = (d_beta_zeta * covariant_theta &
            - d_beta_theta * covariant_zeta) / surface%jacobian

        d_integral_mu0jb_toroidal = field_periods * sum(surface%area_element &
            * d_mu0_toroidal / grad_psi**3) / n_grid
        d_integral_jbsq_toroidal = field_periods * sum(surface%area_element &
            * 2.0_dp * mu0_j_dot_b * d_mu0_toroidal &
            / (b_squared * grad_psi**3)) / n_grid
        d_integral_mu0jb_poloidal = field_periods * sum(surface%area_element &
            * d_mu0_poloidal / grad_psi**3) / n_grid
        d_integral_jbsq_poloidal = field_periods * sum(surface%area_element &
            * 2.0_dp * mu0_j_dot_b * d_mu0_poloidal &
            / (b_squared * grad_psi**3)) / n_grid
    end subroutine flux_beta_integral_derivatives

    subroutine beta_flux_slope_derivative(equilibrium, beta_harmonics, &
            poloidal_flux_slope, toroidal_flux_slope, poloidal_weight, &
            toroidal_weight, theta, zeta, d_beta_theta, d_beta_zeta)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(harmonic_pair_t), intent(in) :: beta_harmonics
        real(dp), intent(in) :: poloidal_flux_slope, toroidal_flux_slope
        real(dp), intent(in) :: poloidal_weight, toroidal_weight
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), allocatable, intent(out) :: d_beta_theta(:, :)
        real(dp), allocatable, intent(out) :: d_beta_zeta(:, :)
        type(harmonic_pair_t) :: d_pair
        real(dp), allocatable :: discard_values(:, :)
        real(dp) :: denominator, d_denominator, scale, mode_m, mode_n
        integer :: idx_m, idx_n, rec_info

        allocate (d_pair%cosine, mold=beta_harmonics%cosine)
        allocate (d_pair%sine, mold=beta_harmonics%sine)
        scale = abs(toroidal_flux_slope) + abs(poloidal_flux_slope)
        do idx_n = 1, size(equilibrium%toroidal_modes)
            mode_n = real(equilibrium%toroidal_modes(idx_n), dp)
            do idx_m = 1, size(equilibrium%poloidal_modes)
                mode_m = real(equilibrium%poloidal_modes(idx_m), dp)
                denominator = two_pi * (mode_m * poloidal_flux_slope &
                    - mode_n * toroidal_flux_slope)
                if (abs(denominator) < 1.0e-10_dp * scale) then
                    d_pair%cosine(1, idx_m, idx_n) = 0.0_dp
                    d_pair%sine(1, idx_m, idx_n) = 0.0_dp
                else
                    d_denominator = two_pi * (poloidal_weight * mode_m &
                        - toroidal_weight * mode_n)
                    d_pair%cosine(1, idx_m, idx_n) = &
                        -beta_harmonics%cosine(1, idx_m, idx_n) &
                        * d_denominator / denominator
                    d_pair%sine(1, idx_m, idx_n) = &
                        -beta_harmonics%sine(1, idx_m, idx_n) &
                        * d_denominator / denominator
                end if
            end do
        end do
        call reconstruct_harmonic_grid(d_pair, 1, &
            equilibrium%poloidal_modes, equilibrium%toroidal_modes, theta, &
            zeta, discard_values, d_beta_theta, d_beta_zeta, rec_info)
    end subroutine beta_flux_slope_derivative

    pure subroutine mercier_flux_slope_gradients(iota_slope, pressure_slope, &
            psi_slope, covariant_zeta, toroidal_flux_slope, &
            current_slope_ratio, volume_curvature, d2v_dpsi2, integral_xi, &
            integral_inverse, integral_bsq, integral_jb, integral_jb_squared, &
            d_integral_mu0jb_toroidal, d_integral_jbsq_toroidal, &
            d_integral_mu0jb_poloidal, d_integral_jbsq_poloidal, &
            d_toroidal_flux_slope, d_poloidal_flux_slope)
        real(dp), intent(in) :: iota_slope, pressure_slope, psi_slope
        real(dp), intent(in) :: covariant_zeta, toroidal_flux_slope
        real(dp), intent(in) :: current_slope_ratio, volume_curvature
        real(dp), intent(in) :: d2v_dpsi2, integral_xi, integral_inverse
        real(dp), intent(in) :: integral_bsq, integral_jb, integral_jb_squared
        real(dp), intent(in) :: d_integral_mu0jb_toroidal
        real(dp), intent(in) :: d_integral_jbsq_toroidal
        real(dp), intent(in) :: d_integral_mu0jb_poloidal
        real(dp), intent(in) :: d_integral_jbsq_poloidal
        real(dp), intent(out) :: d_toroidal_flux_slope, d_poloidal_flux_slope
        real(dp) :: phi, iota_psi, dp_dpsi, current_sign, well_sign
        real(dp) :: well_bracket, d_iota_psi, d_dp_dpsi, d_d2v_dpsi2, d_shear
        real(dp) :: d_xi_t, d_inverse_t, d_bsq_t, d_jb_t, d_jbsq_t
        real(dp) :: d_bracket_t, d_current_t, d_well_t, d_geodesic_t

        phi = toroidal_flux_slope
        iota_psi = iota_slope / psi_slope
        dp_dpsi = mu0 * pressure_slope / psi_slope
        current_sign = -sign(1.0_dp, covariant_zeta)
        well_sign = sign(1.0_dp, psi_slope)
        well_bracket = well_sign * d2v_dpsi2 - dp_dpsi * integral_inverse

        d_iota_psi = -iota_psi / phi
        d_dp_dpsi = -dp_dpsi / phi
        d_shear = 2.0_dp * iota_psi * d_iota_psi &
            / (16.0_dp * acos(-1.0_dp)**2)
        d_d2v_dpsi2 = (volume_curvature / psi_slope**3 &
            - 3.0_dp * d2v_dpsi2 / psi_slope) / two_pi

        d_xi_t = -3.0_dp * integral_xi / phi + d_integral_mu0jb_toroidal &
            + current_slope_ratio / phi * integral_bsq
        d_inverse_t = -integral_inverse / phi
        d_bsq_t = -3.0_dp * integral_bsq / phi
        d_jb_t = -3.0_dp * integral_jb / phi + d_integral_mu0jb_toroidal
        d_jbsq_t = -3.0_dp * integral_jb_squared / phi &
            + d_integral_jbsq_toroidal

        d_current_t = current_sign / two_pi**4 &
            * (d_iota_psi * integral_xi + iota_psi * d_xi_t)
        d_bracket_t = well_sign * d_d2v_dpsi2 &
            - (d_dp_dpsi * integral_inverse + dp_dpsi * d_inverse_t)
        d_well_t = (d_dp_dpsi * well_bracket * integral_bsq &
            + dp_dpsi * d_bracket_t * integral_bsq &
            + dp_dpsi * well_bracket * d_bsq_t) / two_pi**6
        d_geodesic_t = (2.0_dp * integral_jb * d_jb_t &
            - (d_bsq_t * integral_jb_squared + integral_bsq * d_jbsq_t)) &
            / two_pi**6
        d_toroidal_flux_slope = d_shear + d_current_t + d_well_t &
            + d_geodesic_t

        d_poloidal_flux_slope = current_sign / two_pi**4 * iota_psi &
            * d_integral_mu0jb_poloidal &
            + (2.0_dp * integral_jb * d_integral_mu0jb_poloidal &
            - integral_bsq * d_integral_jbsq_poloidal) / two_pi**6
    end subroutine mercier_flux_slope_gradients

    pure function mercier_d_terms(iota_slope, pressure_slope, psi_slope, &
            covariant_zeta, integral_xi, d2v_dpsi2, integral_inverse, &
            integral_bsq, integral_jb, integral_jb_squared) result(d_terms)
        real(dp), intent(in) :: iota_slope, pressure_slope, psi_slope
        real(dp), intent(in) :: covariant_zeta, integral_xi, d2v_dpsi2
        real(dp), intent(in) :: integral_inverse, integral_bsq
        real(dp), intent(in) :: integral_jb, integral_jb_squared
        type(mercier_d_terms_t) :: d_terms
        real(dp) :: iota_psi, dp_dpsi

        iota_psi = iota_slope / psi_slope
        dp_dpsi = mu0 * pressure_slope / psi_slope
        d_terms%shear = iota_psi**2 / (16.0_dp * acos(-1.0_dp)**2)
        d_terms%current = -sign(1.0_dp, covariant_zeta) &
            / two_pi**4 * iota_psi * integral_xi
        d_terms%well = dp_dpsi * (sign(1.0_dp, psi_slope) * d2v_dpsi2 &
            - dp_dpsi * integral_inverse) * integral_bsq / two_pi**6
        d_terms%geodesic = (integral_jb**2 &
            - integral_bsq * integral_jb_squared) / two_pi**6
        d_terms%mercier = d_terms%shear + d_terms%current &
            + d_terms%well + d_terms%geodesic
    end function mercier_d_terms

    pure function mercier_d_terms_gradient(iota_slope, pressure_slope, &
            psi_slope, covariant_zeta, integral_xi, d2v_dpsi2, &
            integral_inverse, integral_bsq) result(grad)
        real(dp), intent(in) :: iota_slope, pressure_slope, psi_slope
        real(dp), intent(in) :: covariant_zeta, integral_xi, d2v_dpsi2
        real(dp), intent(in) :: integral_inverse, integral_bsq
        type(mercier_gradient_t) :: grad
        real(dp) :: dp_dpsi

        dp_dpsi = mu0 * pressure_slope / psi_slope
        grad%d_iota_slope = 2.0_dp * iota_slope &
            / (16.0_dp * acos(-1.0_dp)**2 * psi_slope**2) &
            - sign(1.0_dp, covariant_zeta) * integral_xi &
            / (two_pi**4 * psi_slope)
        grad%d_pressure_slope = (sign(1.0_dp, psi_slope) * d2v_dpsi2 &
            - 2.0_dp * integral_inverse * dp_dpsi) * integral_bsq &
            / two_pi**6 * (mu0 / psi_slope)
    end function mercier_d_terms_gradient

    subroutine load_surface(equilibrium, i, theta, zeta, surface, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: i
        real(dp), intent(in) :: theta(:), zeta(:)
        type(surface_data_t), intent(out) :: surface
        integer, intent(out) :: info

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
        call surface_values(equilibrium%g_st, i, equilibrium, theta, zeta, &
            surface%g_st, info)
        if (info /= mercier_ok) return
        call surface_values(equilibrium%g_sz, i, equilibrium, theta, zeta, &
            surface%g_sz, info)
        if (info /= mercier_ok) return

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

    subroutine solve_beta_derivatives(equilibrium, surface, theta, zeta, &
            covariant_theta_slope, covariant_zeta_slope, pressure_slope, &
            poloidal_flux_slope, toroidal_flux_slope, beta_values, &
            beta_theta, beta_zeta, beta_harmonics)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(surface_data_t), intent(in) :: surface
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(in) :: covariant_theta_slope, covariant_zeta_slope
        real(dp), intent(in) :: pressure_slope, poloidal_flux_slope
        real(dp), intent(in) :: toroidal_flux_slope
        real(dp), allocatable, intent(out) :: beta_values(:, :)
        real(dp), allocatable, intent(out) :: beta_theta(:, :)
        real(dp), allocatable, intent(out) :: beta_zeta(:, :)
        type(harmonic_pair_t), intent(out), optional :: beta_harmonics
        type(harmonic_pair_t) :: beta_pair
        real(dp), allocatable :: rhs(:, :)
        real(dp) :: rhs_cosine(size(equilibrium%poloidal_modes), &
            size(equilibrium%toroidal_modes))
        real(dp) :: rhs_sine(size(equilibrium%poloidal_modes), &
            size(equilibrium%toroidal_modes))
        real(dp) :: denominator, scale
        integer :: mode_m, mode_n, rec_info

        rhs = surface%jacobian * (mu0 * pressure_slope &
            + covariant_zeta_slope * surface%b_zeta &
            + covariant_theta_slope * surface%b_theta)
        call project_harmonic_grid(rhs, equilibrium%poloidal_modes, &
            equilibrium%toroidal_modes, theta, zeta, rhs_cosine, rhs_sine)
        allocate (beta_pair%cosine(1, size(rhs_cosine, 1), &
            size(rhs_cosine, 2)))
        allocate (beta_pair%sine(1, size(rhs_sine, 1), size(rhs_sine, 2)))
        scale = abs(toroidal_flux_slope) + abs(poloidal_flux_slope)
        do mode_n = 1, size(equilibrium%toroidal_modes)
            do mode_m = 1, size(equilibrium%poloidal_modes)
                denominator = two_pi * (real( &
                    equilibrium%poloidal_modes(mode_m), dp) &
                    * poloidal_flux_slope - real( &
                    equilibrium%toroidal_modes(mode_n), dp) &
                    * toroidal_flux_slope)
                if (abs(denominator) < 1.0e-10_dp * scale) then
                    beta_pair%cosine(1, mode_m, mode_n) = 0.0_dp
                    beta_pair%sine(1, mode_m, mode_n) = 0.0_dp
                else
                    beta_pair%sine(1, mode_m, mode_n) = &
                        rhs_cosine(mode_m, mode_n) / denominator
                    beta_pair%cosine(1, mode_m, mode_n) = &
                        -rhs_sine(mode_m, mode_n) / denominator
                end if
            end do
        end do
        call reconstruct_harmonic_grid(beta_pair, 1, &
            equilibrium%poloidal_modes, equilibrium%toroidal_modes, theta, &
            zeta, beta_values, beta_theta, beta_zeta, rec_info)
        if (present(beta_harmonics)) beta_harmonics = beta_pair
    end subroutine solve_beta_derivatives

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

    subroutine filter_gauge_kernel(equilibrium, theta, zeta, &
            poloidal_flux_slope, toroidal_flux_slope, values, filtered)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(in) :: poloidal_flux_slope, toroidal_flux_slope
        real(dp), intent(in) :: values(:, :)
        real(dp), allocatable, intent(out) :: filtered(:, :)
        type(harmonic_pair_t) :: pair
        real(dp), allocatable :: discard_a(:, :), discard_b(:, :)
        real(dp) :: denominator, scale
        integer :: mode_m, mode_n, rec_info

        allocate (pair%cosine(1, size(equilibrium%poloidal_modes), &
            size(equilibrium%toroidal_modes)))
        allocate (pair%sine(1, size(equilibrium%poloidal_modes), &
            size(equilibrium%toroidal_modes)))
        call project_harmonic_grid(values, equilibrium%poloidal_modes, &
            equilibrium%toroidal_modes, theta, zeta, &
            pair%cosine(1, :, :), pair%sine(1, :, :))
        scale = abs(toroidal_flux_slope) + abs(poloidal_flux_slope)
        do mode_n = 1, size(equilibrium%toroidal_modes)
            do mode_m = 1, size(equilibrium%poloidal_modes)
                denominator = two_pi * (real( &
                    equilibrium%poloidal_modes(mode_m), dp) &
                    * poloidal_flux_slope - real( &
                    equilibrium%toroidal_modes(mode_n), dp) &
                    * toroidal_flux_slope)
                if (abs(denominator) < 1.0e-10_dp * scale) then
                    pair%cosine(1, mode_m, mode_n) = 0.0_dp
                    pair%sine(1, mode_m, mode_n) = 0.0_dp
                end if
            end do
        end do
        call reconstruct_harmonic_grid(pair, 1, &
            equilibrium%poloidal_modes, equilibrium%toroidal_modes, &
            theta, zeta, filtered, discard_a, discard_b, rec_info)
    end subroutine filter_gauge_kernel

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

    subroutine beta_from_positions(equilibrium, xhat_s, yhat_s, zhat_s, &
            i, theta, zeta, surface, beta_positions, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        type(harmonic_pair_t), intent(in) :: xhat_s, yhat_s, zhat_s
        integer, intent(in) :: i
        real(dp), intent(in) :: theta(:), zeta(:)
        type(surface_data_t), intent(in) :: surface
        real(dp), allocatable, intent(out) :: beta_positions(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: es_x(:, :), es_y(:, :), es_z(:, :)
        real(dp), allocatable :: eu_x(:, :), eu_y(:, :), eu_z(:, :)
        real(dp), allocatable :: ev_x(:, :), ev_y(:, :), ev_z(:, :)
        real(dp), allocatable :: discard_a(:, :), discard_b(:, :)

        call surface_values(xhat_s, i, equilibrium, theta, zeta, es_x, &
            info)
        if (info /= mercier_ok) return
        call surface_values(yhat_s, i, equilibrium, theta, zeta, es_y, &
            info)
        if (info /= mercier_ok) return
        call surface_values(zhat_s, i, equilibrium, theta, zeta, es_z, &
            info)
        if (info /= mercier_ok) return
        call surface_derivatives(equilibrium%xhat, i, equilibrium, theta, &
            zeta, discard_a, eu_x, ev_x, info)
        if (info /= mercier_ok) return
        call surface_derivatives(equilibrium%yhat, i, equilibrium, theta, &
            zeta, discard_a, eu_y, ev_y, info)
        if (info /= mercier_ok) return
        call surface_derivatives(equilibrium%zhat, i, equilibrium, theta, &
            zeta, discard_b, eu_z, ev_z, info)
        if (info /= mercier_ok) return
        beta_positions = (es_x * eu_x + es_y * eu_y + es_z * eu_z) &
            * surface%b_theta &
            + (es_x * ev_x + es_y * ev_y + es_z * ev_z) * surface%b_zeta
    end subroutine beta_from_positions

    subroutine surface_derivatives(pair, i, equilibrium, theta, zeta, &
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
    end subroutine surface_derivatives

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
        allocate (result%d_mercier_d_iota_slope(ns))
        allocate (result%d_mercier_d_pressure_slope(ns))
        allocate (result%d_mercier_d_pressure_slope_full(ns))
        allocate (result%d_mercier_d_toroidal_flux_slope(ns))
        allocate (result%d_mercier_d_poloidal_flux_slope(ns))
        allocate (result%boozer_deviation(ns))
        allocate (result%force_balance_residual(ns))
        allocate (result%jacobian_identity_deviation(ns))
        allocate (result%beta_chart_deviation(ns))
    end subroutine allocate_result

end module mercier_diagnostic
