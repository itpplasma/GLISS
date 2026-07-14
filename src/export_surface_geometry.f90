module export_surface_geometry
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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

    real(dp), parameter, public :: two_pi = 2.0_dp * acos(-1.0_dp)
    real(dp), parameter, public :: mu0 = 2.0_dp * two_pi * 1.0e-7_dp

    type, public :: surface_data_t
        real(dp), allocatable :: jacobian(:, :), g_tt(:, :), g_tz(:, :)
        real(dp), allocatable :: g_zz(:, :), b_theta(:, :), b_zeta(:, :)
        real(dp), allocatable :: g_st(:, :), g_sz(:, :)
        real(dp), allocatable :: mod_b(:, :)
        real(dp), allocatable :: area_element(:, :)
    end type surface_data_t

    type, public :: surface_profiles_t
        real(dp) :: flux_slope, poloidal_slope
        real(dp) :: flux_curvature, poloidal_curvature
        real(dp) :: covariant_theta, covariant_zeta
        real(dp) :: covariant_theta_slope, covariant_zeta_slope
        real(dp) :: pressure_slope
    end type surface_profiles_t

    public :: build_angular_grids
    public :: build_kernel_geometry
    public :: beta_from_positions
    public :: differentiate_pair
    public :: grid_mean
    public :: load_surface
    public :: solve_beta_derivatives
    public :: surface_derivatives
    public :: surface_values

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
        real(dp), allocatable :: grad_s2(:, :)
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
        grad_s2 = (surface%g_tt * surface%g_zz - surface%g_tz**2) &
            / surface%jacobian**2
        if (.not. all(ieee_is_finite(grad_s2)) &
            .or. any(grad_s2 <= 0.0_dp)) then
            info = mercier_invalid_input
            return
        end if
        fields(:, :, 9) = grad_s2
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
        grad_s2 = (surface%g_tt * surface%g_zz - surface%g_tz**2) &
            / surface%jacobian**2
        if (.not. all(ieee_is_finite(grad_s2)) &
            .or. any(grad_s2 <= 0.0_dp)) then
            info = mercier_invalid_input
            return
        end if
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

end module export_surface_geometry
