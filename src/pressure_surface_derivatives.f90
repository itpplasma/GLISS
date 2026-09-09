! Contract v1: pressure samples in Pa on fixed normalized-flux nodes s in [0,1].
! Geometry, flux profiles, harmonic topology, angular grid, gamma and coordinate
! are fixed. This differentiates the production surface fields/drive and gamma*p;
! it does not differentiate force balance or assemble global K/M derivatives.
module pressure_surface_derivatives
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use export_surface_geometry, only: build_surface_kernel_fields, &
        mercier_ok, mu0, solve_beta_derivatives_modes, surface_data_t, &
        surface_profiles_t
    use gvec_cas3d_reconstruction, only: project_harmonic_grid, &
        reconstruct_harmonic_grid, reconstruction_ok
    use gvec_cas3d_types, only: harmonic_pair_t
    use radial_cubic_spline, only: build_radial_cubic_spline_grid, &
        evaluate_radial_cubic_spline_field, fit_radial_cubic_spline_field, &
        radial_cubic_spline_field_t, radial_cubic_spline_grid_t, &
        radial_cubic_spline_ok
    implicit none
    private

    integer, parameter, public :: pressure_derivative_ok = 0
    integer, parameter, public :: pressure_derivative_invalid = -1
    integer, parameter, public :: pressure_derivative_allocation = -2
    integer, parameter, public :: pressure_derivative_contract_version = 1

    ! fields_slope and drive_slope are derivatives with respect to p_s, not s.
    ! gamma_pressure_weights and slope_weights act on the input pressure samples.
    ! Storage for the cardinal spline solve is quadratic in sample count.
    type, public :: pressure_surface_response_t
        logical :: ready = .false.
        real(dp) :: gamma_pressure = 0.0_dp
        real(dp), allocatable :: fields(:, :, :), drive(:, :)
        real(dp), allocatable :: fields_slope(:, :, :), drive_slope(:, :)
        real(dp), allocatable :: gamma_pressure_weights(:), slope_weights(:)
    end type pressure_surface_response_t

    public :: build_pressure_surface_response
    public :: pressure_surface_jvp, pressure_surface_vjp

contains

    subroutine build_pressure_surface_response(m, n, chart, surface, profiles, &
            jacobian_slope, theta, zeta, nodes, pressure, coordinate, gamma, &
            response, info)
        integer, intent(in) :: m(:), n(:)
        logical, intent(in) :: chart
        type(surface_data_t), intent(in) :: surface
        type(surface_profiles_t), intent(in) :: profiles
        real(dp), intent(in) :: jacobian_slope(:, :), theta(:), zeta(:)
        real(dp), intent(in) :: nodes(:), pressure(:), coordinate, gamma
        type(pressure_surface_response_t), intent(out) :: response
        integer, intent(out) :: info
        type(surface_profiles_t) :: active
        real(dp) :: value, slope
        integer :: status

        info = pressure_derivative_invalid
        if (.not. ieee_is_finite(gamma)) return
        if (gamma <= 0.0_dp) return
        call pressure_sample_map(nodes, pressure, coordinate, value, slope, &
            response%gamma_pressure_weights, response%slope_weights, status)
        if (status /= pressure_derivative_ok) then
            info = status
            return
        end if
        ! An open positive-pressure domain is required for material derivatives.
        if (value <= 0.0_dp) return
        active = profiles
        active%pressure_slope = slope
        allocate (response%fields(size(theta), size(zeta), 13), &
            response%drive(size(theta), size(zeta)), &
            response%fields_slope(size(theta), size(zeta), 13), &
            response%drive_slope(size(theta), size(zeta)), stat=status)
        if (status /= 0) then
            info = pressure_derivative_allocation
            return
        end if
        call build_surface_kernel_fields(m, n, chart, surface, active, &
            jacobian_slope, theta, zeta, response%fields, response%drive, status)
        if (status /= mercier_ok) return
        call surface_slope_response(m, n, chart, surface, active, &
            jacobian_slope, theta, zeta, response%fields, &
            response%fields_slope, response%drive_slope, status)
        if (status /= pressure_derivative_ok) return
        response%gamma_pressure = gamma * value
        response%gamma_pressure_weights = gamma * response%gamma_pressure_weights
        if (.not. ieee_is_finite(response%gamma_pressure)) return
        if (.not. all(ieee_is_finite(response%fields))) return
        if (.not. all(ieee_is_finite(response%drive))) return
        if (.not. all(ieee_is_finite(response%gamma_pressure_weights))) return
        response%ready = .true.
        info = pressure_derivative_ok
    end subroutine build_pressure_surface_response

    subroutine pressure_sample_map(nodes, pressure, coordinate, value, slope, &
            value_weights, slope_weights, info)
        real(dp), intent(in) :: nodes(:), pressure(:), coordinate
        real(dp), intent(out) :: value, slope
        real(dp), allocatable, intent(out) :: value_weights(:), slope_weights(:)
        integer, intent(out) :: info
        type(radial_cubic_spline_grid_t) :: grid
        type(radial_cubic_spline_field_t) :: spline
        real(dp), allocatable :: samples(:, :), values(:), slopes(:), second(:)
        integer :: count, column, status

        info = pressure_derivative_invalid
        value = 0.0_dp
        slope = 0.0_dp
        count = size(nodes)
        if (size(pressure) /= count) return
        if (.not. all(ieee_is_finite(pressure))) return
        if (any(pressure <= 0.0_dp)) return
        call build_radial_cubic_spline_grid(nodes, 0.0_dp, 1.0_dp, grid, status)
        if (status /= radial_cubic_spline_ok) return
        ! Cardinal data produce the exact linear operator of the SAME spline
        ! solver used by primitive_equilibrium_spline, not a perturbed solve.
        allocate (samples(count, count + 1), values(count + 1), &
            slopes(count + 1), second(count + 1), stat=status)
        if (status /= 0) then
            info = pressure_derivative_allocation
            return
        end if
        samples = 0.0_dp
        samples(:, 1) = pressure
        do column = 1, count
            samples(column, column + 1) = 1.0_dp
        end do
        call fit_radial_cubic_spline_field(grid, samples, spline, status)
        if (status /= radial_cubic_spline_ok) return
        call evaluate_radial_cubic_spline_field(grid, spline, coordinate, &
            values, slopes, second, status)
        if (status /= radial_cubic_spline_ok) return
        value = values(1)
        slope = slopes(1)
        value_weights = values(2:)
        slope_weights = slopes(2:)
        info = pressure_derivative_ok
    end subroutine pressure_sample_map

    subroutine surface_slope_response(m, n, chart, surface, profiles, &
            jacobian_slope, theta, zeta, fields, fields_slope, drive_slope, info)
        integer, intent(in) :: m(:), n(:)
        logical, intent(in) :: chart
        type(surface_data_t), intent(in) :: surface
        type(surface_profiles_t), intent(in) :: profiles
        real(dp), intent(in) :: jacobian_slope(:, :), theta(:), zeta(:)
        real(dp), intent(in) :: fields(:, :, :)
        real(dp), intent(out) :: fields_slope(:, :, :), drive_slope(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: beta(:, :), beta_theta(:, :), beta_zeta(:, :)
        integer :: status

        info = pressure_derivative_invalid
        ! Fixed denominator: the beta derivative is the production inverse MDE
        ! applied to the pressure-only forcing. This also rejects a compatible
        ! primal whose pressure tangent violates a resonant compatibility rule.
        call solve_beta_derivatives_modes(m, n, surface, theta, zeta, &
            0.0_dp, 0.0_dp, 1.0_dp, profiles%poloidal_slope, &
            profiles%flux_slope, beta, beta_theta, beta_zeta, info=status)
        if (status /= mercier_ok) return
        fields_slope = 0.0_dp
        fields_slope(:, :, 10) = (beta_zeta * fields(:, :, 6) &
            - beta_theta * fields(:, :, 5)) / surface%jacobian
        fields_slope(:, :, 11) = mu0
        if (.not. chart) fields_slope(:, :, 13) = beta
        drive_slope = 2.0_dp * (fields(:, :, 10) * fields_slope(:, :, 10) &
            + mu0**2 * profiles%pressure_slope * fields(:, :, 9)) &
            / (surface%mod_b**2 * fields(:, :, 9)) &
            + (profiles%flux_curvature * beta_zeta &
            + profiles%poloidal_curvature * beta_theta) / surface%jacobian &
            - mu0 * jacobian_slope / surface%jacobian
        if (chart) then
            call chart_slope_response(m, n, surface, profiles, theta, zeta, &
                beta_theta, beta_zeta, drive_slope, status)
            if (status /= pressure_derivative_ok) return
        end if
        if (.not. all(ieee_is_finite(fields_slope))) return
        if (.not. all(ieee_is_finite(drive_slope))) return
        info = pressure_derivative_ok
    end subroutine surface_slope_response

    subroutine chart_slope_response(m, n, surface, profiles, theta, zeta, &
            beta_theta, beta_zeta, drive_slope, info)
        integer, intent(in) :: m(:), n(:)
        type(surface_data_t), intent(in) :: surface
        type(surface_profiles_t), intent(in) :: profiles
        real(dp), intent(in) :: theta(:), zeta(:), beta_theta(:, :), beta_zeta(:, :)
        real(dp), intent(inout) :: drive_slope(:, :)
        integer, intent(out) :: info
        type(harmonic_pair_t) :: pair
        real(dp), allocatable :: operand(:, :), values(:, :), dt(:, :), dz(:, :)
        real(dp) :: cosine(size(m), size(n)), sine(size(m), size(n))
        integer :: status

        info = pressure_derivative_invalid
        operand = -(beta_theta * (surface%g_sz * surface%g_tz &
            - surface%g_st * surface%g_zz) &
            + beta_zeta * (surface%g_st * surface%g_tz &
            - surface%g_sz * surface%g_tt)) &
            / (surface%g_tt * surface%g_zz - surface%g_tz**2)
        call project_harmonic_grid(operand, m, n, theta, zeta, cosine, sine)
        allocate (pair%cosine(1, size(m), size(n)), &
            pair%sine(1, size(m), size(n)), stat=status)
        if (status /= 0) return
        pair%cosine(1, :, :) = cosine
        pair%sine(1, :, :) = sine
        call reconstruct_harmonic_grid(pair, 1, m, n, theta, zeta, &
            values, dt, dz, status)
        if (status /= reconstruction_ok) return
        drive_slope = drive_slope + (profiles%poloidal_slope * dt &
            + profiles%flux_slope * dz) / surface%jacobian
        info = pressure_derivative_ok
    end subroutine chart_slope_response

    subroutine pressure_surface_jvp(response, tangent, fields, drive, &
            gamma_pressure, info)
        type(pressure_surface_response_t), intent(in) :: response
        real(dp), intent(in) :: tangent(:)
        real(dp), allocatable, intent(out) :: fields(:, :, :), drive(:, :)
        real(dp), intent(out) :: gamma_pressure
        integer, intent(out) :: info
        real(dp) :: slope

        info = pressure_derivative_invalid
        gamma_pressure = 0.0_dp
        if (.not. response%ready) return
        if (size(tangent) /= size(response%slope_weights)) return
        if (.not. all(ieee_is_finite(tangent))) return
        slope = dot_product(response%slope_weights, tangent)
        gamma_pressure = dot_product(response%gamma_pressure_weights, tangent)
        fields = slope * response%fields_slope
        drive = slope * response%drive_slope
        if (.not. ieee_is_finite(gamma_pressure) &
            .or. .not. all(ieee_is_finite(fields)) &
            .or. .not. all(ieee_is_finite(drive))) then
            deallocate (fields, drive)
            gamma_pressure = 0.0_dp
            return
        end if
        info = pressure_derivative_ok
    end subroutine pressure_surface_jvp

    subroutine pressure_surface_vjp(response, fields, drive, gamma_pressure, &
            gradient, info)
        type(pressure_surface_response_t), intent(in) :: response
        real(dp), intent(in) :: fields(:, :, :), drive(:, :), gamma_pressure
        real(dp), allocatable, intent(out) :: gradient(:)
        integer, intent(out) :: info
        real(dp) :: weight

        info = pressure_derivative_invalid
        if (.not. response%ready) return
        if (any(shape(fields) /= shape(response%fields_slope))) return
        if (any(shape(drive) /= shape(response%drive_slope))) return
        if (.not. all(ieee_is_finite(fields))) return
        if (.not. all(ieee_is_finite(drive))) return
        if (.not. ieee_is_finite(gamma_pressure)) return
        weight = sum(fields * response%fields_slope) &
            + sum(drive * response%drive_slope)
        gradient = weight * response%slope_weights &
            + gamma_pressure * response%gamma_pressure_weights
        if (.not. all(ieee_is_finite(gradient))) then
            deallocate (gradient)
            return
        end if
        info = pressure_derivative_ok
    end subroutine pressure_surface_vjp
end module pressure_surface_derivatives
