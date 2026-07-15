module field_profile_identities
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gvec_cas3d_reconstruction, only: reconstruct_harmonic_grid, &
        reconstruction_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t, &
        harmonic_pair_t, radial_grid_half
    use nonuniform_derivative, only: first_derivative_nonuniform
    implicit none
    private

    integer, parameter, public :: field_profile_identities_ok = 0
    integer, parameter, public :: field_profile_identities_invalid_input = 1
    integer, parameter, public :: field_profile_identities_reconstruction_error = 2

    real(dp), parameter :: mu0 = 4.0e-7_dp * acos(-1.0_dp)

    type, public :: field_profile_identity_result_t
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: field_toroidal_flux_slope(:)
        real(dp), allocatable :: field_poloidal_flux_slope(:)
        real(dp), allocatable :: field_covariant_theta(:)
        real(dp), allocatable :: field_covariant_zeta(:)
        real(dp), allocatable :: toroidal_flux_deviation(:)
        real(dp), allocatable :: poloidal_flux_deviation(:)
        real(dp), allocatable :: covariant_theta_deviation(:)
        real(dp), allocatable :: covariant_zeta_deviation(:)
        real(dp), allocatable :: iota_flux_deviation(:)
        real(dp), allocatable :: ampere_theta_deviation(:)
        real(dp), allocatable :: ampere_zeta_deviation(:)
        real(dp), allocatable :: exported_jacobian_deviation(:)
        real(dp), allocatable :: general_force_balance_deviation(:)
    end type field_profile_identity_result_t

    public :: compute_field_profile_identities

contains

    pure subroutine compute_field_profile_identities(equilibrium, n_theta, &
            n_zeta, result, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: n_theta, n_zeta
        type(field_profile_identity_result_t), intent(out) :: result
        integer, intent(out) :: info
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp), allocatable :: jacobian(:, :, :)
        real(dp), allocatable :: mod_b_squared(:, :, :)
        real(dp), allocatable :: contra_theta(:, :, :), contra_zeta(:, :, :)
        real(dp), allocatable :: b_theta(:, :, :), b_zeta(:, :, :)
        real(dp), allocatable :: beta_theta(:, :, :), beta_zeta(:, :, :)
        real(dp), allocatable :: flux_toroidal_slope(:)
        real(dp), allocatable :: flux_poloidal_slope(:)
        real(dp), allocatable :: b_theta_slope(:), b_zeta_slope(:)
        real(dp), allocatable :: pressure_slope(:)
        real(dp), allocatable :: field_theta_slope(:, :, :)
        real(dp), allocatable :: field_zeta_slope(:, :, :)
        real(dp) :: periods
        integer :: ns, i, j, k

        info = field_profile_identities_invalid_input
        if (n_theta < 8 .or. n_zeta < 8) return
        if (equilibrium%radial_grid /= radial_grid_half) return
        ns = size(equilibrium%s)
        if (ns < 3) return
        call allocate_result(ns, result)
        call build_angular_grids(n_theta, n_zeta, theta, zeta)
        allocate (jacobian(n_theta, n_zeta, ns))
        allocate (mod_b_squared(n_theta, n_zeta, ns))
        allocate (contra_theta(n_theta, n_zeta, ns))
        allocate (contra_zeta(n_theta, n_zeta, ns))
        allocate (b_theta(n_theta, n_zeta, ns))
        allocate (b_zeta(n_theta, n_zeta, ns))
        allocate (beta_theta(n_theta, n_zeta, ns))
        allocate (beta_zeta(n_theta, n_zeta, ns))
        do i = 1, ns
            call reconstruct_surface(equilibrium, i, theta, zeta, &
                jacobian(:, :, i), mod_b_squared(:, :, i), &
                contra_theta(:, :, i), contra_zeta(:, :, i), &
                b_theta(:, :, i), b_zeta(:, :, i), beta_theta(:, :, i), &
                beta_zeta(:, :, i), info)
            if (info /= field_profile_identities_ok) return
        end do

        allocate (flux_toroidal_slope(ns), flux_poloidal_slope(ns))
        allocate (b_theta_slope(ns), b_zeta_slope(ns), pressure_slope(ns))
        call first_derivative_nonuniform(equilibrium%s, &
            equilibrium%toroidal_flux, flux_toroidal_slope)
        call first_derivative_nonuniform(equilibrium%s, &
            equilibrium%poloidal_flux, flux_poloidal_slope)
        call first_derivative_nonuniform(equilibrium%s, &
            equilibrium%b_theta_average, b_theta_slope)
        call first_derivative_nonuniform(equilibrium%s, &
            equilibrium%b_zeta_average, b_zeta_slope)
        call first_derivative_nonuniform(equilibrium%s, &
            equilibrium%pressure, pressure_slope)
        allocate (field_theta_slope(n_theta, n_zeta, ns))
        allocate (field_zeta_slope(n_theta, n_zeta, ns))
        do k = 1, n_zeta
            do j = 1, n_theta
                call first_derivative_nonuniform(equilibrium%s, &
                    b_theta(j, k, :), field_theta_slope(j, k, :))
                call first_derivative_nonuniform(equilibrium%s, &
                    b_zeta(j, k, :), field_zeta_slope(j, k, :))
            end do
        end do

        periods = real(equilibrium%field_periods, dp)
        do i = 1, ns
            call summarize_surface(equilibrium, i, periods, jacobian(:, :, i), &
                mod_b_squared(:, :, i), contra_theta(:, :, i), &
                contra_zeta(:, :, i), b_theta(:, :, i), b_zeta(:, :, i), &
                beta_theta(:, :, i), beta_zeta(:, :, i), &
                field_theta_slope(:, :, i), field_zeta_slope(:, :, i), &
                flux_toroidal_slope, flux_poloidal_slope, b_theta_slope, &
                b_zeta_slope, pressure_slope, result)
        end do
        result%s = equilibrium%s
        info = field_profile_identities_ok
    end subroutine compute_field_profile_identities

    pure subroutine reconstruct_surface(equilibrium, radial_surface, theta, &
            zeta, jacobian, mod_b_squared, contra_theta, contra_zeta, &
            b_theta, b_zeta, beta_theta, beta_zeta, info)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: radial_surface
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), intent(out) :: jacobian(:, :), mod_b_squared(:, :)
        real(dp), intent(out) :: contra_theta(:, :), contra_zeta(:, :)
        real(dp), intent(out) :: b_theta(:, :), b_zeta(:, :)
        real(dp), intent(out) :: beta_theta(:, :), beta_zeta(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: g_tt(:, :), g_tz(:, :), g_zz(:, :)
        real(dp), allocatable :: g_st(:, :), g_sz(:, :)
        real(dp), allocatable :: mod_b(:, :)
        real(dp), allocatable :: jacobian_values(:, :)
        real(dp), allocatable :: contra_theta_values(:, :)
        real(dp), allocatable :: contra_zeta_values(:, :)
        real(dp), allocatable :: g_st_theta(:, :), g_st_zeta(:, :)
        real(dp), allocatable :: g_sz_theta(:, :), g_sz_zeta(:, :)
        real(dp), allocatable :: contra_theta_theta(:, :)
        real(dp), allocatable :: contra_theta_zeta(:, :)
        real(dp), allocatable :: contra_zeta_theta(:, :)
        real(dp), allocatable :: contra_zeta_zeta(:, :)

        call reconstruct_values(equilibrium%jacobian, radial_surface, &
            equilibrium, theta, zeta, jacobian_values, info)
        if (info /= field_profile_identities_ok) return
        call reconstruct_values(equilibrium%mod_b, radial_surface, &
            equilibrium, theta, zeta, mod_b, info)
        if (info /= field_profile_identities_ok) return
        call reconstruct_values(equilibrium%g_tt, radial_surface, &
            equilibrium, theta, zeta, g_tt, info)
        if (info /= field_profile_identities_ok) return
        call reconstruct_values(equilibrium%g_tz, radial_surface, &
            equilibrium, theta, zeta, g_tz, info)
        if (info /= field_profile_identities_ok) return
        call reconstruct_values(equilibrium%g_zz, radial_surface, &
            equilibrium, theta, zeta, g_zz, info)
        if (info /= field_profile_identities_ok) return
        call reconstruct_with_derivatives(equilibrium%g_st, radial_surface, &
            equilibrium, theta, zeta, g_st, g_st_theta, g_st_zeta, info)
        if (info /= field_profile_identities_ok) return
        call reconstruct_with_derivatives(equilibrium%g_sz, radial_surface, &
            equilibrium, theta, zeta, g_sz, g_sz_theta, g_sz_zeta, info)
        if (info /= field_profile_identities_ok) return
        call reconstruct_with_derivatives( &
            equilibrium%b_contravariant_theta, radial_surface, equilibrium, &
            theta, zeta, contra_theta_values, contra_theta_theta, &
            contra_theta_zeta, info)
        if (info /= field_profile_identities_ok) return
        call reconstruct_with_derivatives( &
            equilibrium%b_contravariant_zeta, radial_surface, equilibrium, &
            theta, zeta, contra_zeta_values, contra_zeta_theta, &
            contra_zeta_zeta, info)
        if (info /= field_profile_identities_ok) return
        jacobian = jacobian_values
        contra_theta = contra_theta_values
        contra_zeta = contra_zeta_values
        b_theta = g_tt * contra_theta + g_tz * contra_zeta
        b_zeta = g_tz * contra_theta + g_zz * contra_zeta
        mod_b_squared = mod_b**2
        beta_theta = g_st_theta * contra_theta &
            + g_st * contra_theta_theta + g_sz_theta * contra_zeta &
            + g_sz * contra_zeta_theta
        beta_zeta = g_st_zeta * contra_theta &
            + g_st * contra_theta_zeta + g_sz_zeta * contra_zeta &
            + g_sz * contra_zeta_zeta
    end subroutine reconstruct_surface

    pure subroutine reconstruct_values(pair, radial_surface, equilibrium, &
            theta, zeta, values, info)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: radial_surface
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), allocatable, intent(out) :: values(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: derivative_theta(:, :)
        real(dp), allocatable :: derivative_zeta(:, :)

        call reconstruct_harmonic_grid(pair, radial_surface, &
            equilibrium%poloidal_modes, equilibrium%toroidal_modes, theta, &
            zeta, values, derivative_theta, derivative_zeta, info)
        if (info /= reconstruction_ok) then
            info = field_profile_identities_reconstruction_error
        else
            info = field_profile_identities_ok
        end if
    end subroutine reconstruct_values

    pure subroutine reconstruct_with_derivatives(pair, radial_surface, &
            equilibrium, theta, zeta, values, derivative_theta, &
            derivative_zeta, info)
        type(harmonic_pair_t), intent(in) :: pair
        integer, intent(in) :: radial_surface
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: theta(:), zeta(:)
        real(dp), allocatable, intent(out) :: values(:, :)
        real(dp), allocatable, intent(out) :: derivative_theta(:, :)
        real(dp), allocatable, intent(out) :: derivative_zeta(:, :)
        integer, intent(out) :: info

        call reconstruct_harmonic_grid(pair, radial_surface, &
            equilibrium%poloidal_modes, equilibrium%toroidal_modes, theta, &
            zeta, values, derivative_theta, derivative_zeta, info)
        if (info /= reconstruction_ok) then
            info = field_profile_identities_reconstruction_error
        else
            info = field_profile_identities_ok
        end if
    end subroutine reconstruct_with_derivatives

    pure subroutine summarize_surface(equilibrium, i, periods, jacobian, &
            mod_b_squared, contra_theta, contra_zeta, b_theta, b_zeta, &
            beta_theta, beta_zeta, field_theta_slope, &
            field_zeta_slope, flux_toroidal_slope, flux_poloidal_slope, &
            b_theta_slope, b_zeta_slope, pressure_slope, result)
        type(gvec_cas3d_equilibrium_t), intent(in) :: equilibrium
        integer, intent(in) :: i
        real(dp), intent(in) :: periods
        real(dp), intent(in) :: jacobian(:, :), b_theta(:, :), b_zeta(:, :)
        real(dp), intent(in) :: mod_b_squared(:, :)
        real(dp), intent(in) :: contra_theta(:, :), contra_zeta(:, :)
        real(dp), intent(in) :: beta_theta(:, :), beta_zeta(:, :)
        real(dp), intent(in) :: field_theta_slope(:, :)
        real(dp), intent(in) :: field_zeta_slope(:, :)
        real(dp), intent(in) :: flux_toroidal_slope(:)
        real(dp), intent(in) :: flux_poloidal_slope(:)
        real(dp), intent(in) :: b_theta_slope(:), b_zeta_slope(:)
        real(dp), intent(in) :: pressure_slope(:)
        type(field_profile_identity_result_t), intent(inout) :: result
        real(dp), allocatable :: residual(:, :), terms(:, :)
        real(dp) :: covariant_scale, current_scale, flux_scale
        real(dp) :: poloidal_field_slope, toroidal_field_slope
        integer :: column, row

        toroidal_field_slope = -flux_toroidal_slope(i)
        poloidal_field_slope = -flux_poloidal_slope(i) / periods
        flux_scale = max(maxval(abs(flux_toroidal_slope)), &
            maxval(abs(flux_poloidal_slope)), tiny(1.0_dp))
        covariant_scale = max(maxval(abs(equilibrium%b_theta_average)), &
            maxval(abs(equilibrium%b_zeta_average)), tiny(1.0_dp))
        current_scale = max(maxval(abs(b_theta_slope)), &
            maxval(abs(b_zeta_slope)), tiny(1.0_dp))
        allocate (residual(size(jacobian, 1), size(jacobian, 2)), &
            terms(size(jacobian, 1), size(jacobian, 2)))
        result%field_toroidal_flux_slope(i) = &
            -grid_product_mean(jacobian, contra_zeta)
        result%field_poloidal_flux_slope(i) = &
            -periods * grid_product_mean(jacobian, contra_theta)
        result%field_covariant_theta(i) = grid_mean(b_theta)
        result%field_covariant_zeta(i) = grid_mean(b_zeta)
        result%toroidal_flux_deviation(i) = &
            relative_product_deviation(jacobian, contra_zeta, -1.0_dp, &
            flux_toroidal_slope(i), flux_scale)
        result%poloidal_flux_deviation(i) = &
            relative_product_deviation(jacobian, contra_theta, -periods, &
            flux_poloidal_slope(i), flux_scale)
        result%covariant_theta_deviation(i) = relative_grid_deviation( &
            b_theta, equilibrium%b_theta_average(i), covariant_scale)
        result%covariant_zeta_deviation(i) = relative_grid_deviation( &
            b_zeta, equilibrium%b_zeta_average(i), covariant_scale)
        result%iota_flux_deviation(i) = relative_scalar_deviation( &
            flux_poloidal_slope(i), equilibrium%rotational_transform(i) &
            * flux_toroidal_slope(i), flux_scale)
        result%ampere_theta_deviation(i) = relative_scalar_deviation( &
            grid_difference_mean(beta_zeta, field_zeta_slope), &
            -b_zeta_slope(i), &
            current_scale)
        result%ampere_zeta_deviation(i) = relative_scalar_deviation( &
            grid_difference_mean(field_theta_slope, beta_theta), &
            b_theta_slope(i), &
            current_scale)
        do column = 1, size(jacobian, 2)
            do row = 1, size(jacobian, 1)
                residual(row, column) = mod_b_squared(row, column) &
                    * jacobian(row, column) - toroidal_field_slope &
                    * equilibrium%b_zeta_average(i) - poloidal_field_slope &
                    * equilibrium%b_theta_average(i)
                terms(row, column) = abs(mod_b_squared(row, column) &
                    * jacobian(row, column)) + abs(toroidal_field_slope &
                    * equilibrium%b_zeta_average(i)) &
                    + abs(poloidal_field_slope &
                    * equilibrium%b_theta_average(i))
            end do
        end do
        result%exported_jacobian_deviation(i) = scaled_grid_residual( &
            residual, terms)
        do column = 1, size(jacobian, 2)
            do row = 1, size(jacobian, 1)
                residual(row, column) = mu0 * pressure_slope(i) &
                    * jacobian(row, column) + toroidal_field_slope &
                    * b_zeta_slope(i) + poloidal_field_slope &
                    * b_theta_slope(i) - toroidal_field_slope &
                    * beta_zeta(row, column) - poloidal_field_slope &
                    * beta_theta(row, column)
                terms(row, column) = abs(mu0 * pressure_slope(i) &
                    * jacobian(row, column)) + abs(toroidal_field_slope &
                    * b_zeta_slope(i)) + abs(poloidal_field_slope &
                    * b_theta_slope(i)) + abs(toroidal_field_slope &
                    * beta_zeta(row, column)) + abs(poloidal_field_slope &
                    * beta_theta(row, column))
            end do
        end do
        result%general_force_balance_deviation(i) = scaled_grid_residual( &
            residual, terms)
    end subroutine summarize_surface

    pure function relative_grid_deviation(values, reference, scale) &
            result(value)
        real(dp), intent(in) :: values(:, :), reference, scale
        real(dp) :: value

        value = maxval(abs(values - reference)) / scale
    end function relative_grid_deviation

    pure function relative_scalar_deviation(value_a, value_b, scale) &
            result(value)
        real(dp), intent(in) :: value_a, value_b, scale
        real(dp) :: value

        value = abs(value_a - value_b) / scale
    end function relative_scalar_deviation

    pure function scaled_grid_residual(residual, terms) result(value)
        real(dp), intent(in) :: residual(:, :), terms(:, :)
        real(dp) :: value

        value = maxval(abs(residual)) &
            / max(maxval(terms), tiny(1.0_dp))
    end function scaled_grid_residual

    pure function grid_mean(values) result(value)
        real(dp), intent(in) :: values(:, :)
        real(dp) :: value

        value = sum(values) / real(size(values), dp)
    end function grid_mean

    pure function grid_product_mean(first, second) result(value)
        real(dp), intent(in) :: first(:, :), second(:, :)
        real(dp) :: value
        integer :: column, row

        value = 0.0_dp
        do column = 1, size(first, 2)
            do row = 1, size(first, 1)
                value = value + first(row, column) * second(row, column)
            end do
        end do
        value = value / real(size(first), dp)
    end function grid_product_mean

    pure function grid_difference_mean(first, second) result(value)
        real(dp), intent(in) :: first(:, :), second(:, :)
        real(dp) :: value
        integer :: column, row

        value = 0.0_dp
        do column = 1, size(first, 2)
            do row = 1, size(first, 1)
                value = value + first(row, column) - second(row, column)
            end do
        end do
        value = value / real(size(first), dp)
    end function grid_difference_mean

    pure function relative_product_deviation(first, second, multiplier, &
            reference, scale) result(value)
        real(dp), intent(in) :: first(:, :), second(:, :)
        real(dp), intent(in) :: multiplier, reference, scale
        real(dp) :: value
        integer :: column, row

        value = 0.0_dp
        do column = 1, size(first, 2)
            do row = 1, size(first, 1)
                value = max(value, abs(multiplier * first(row, column) &
                    * second(row, column) - reference))
            end do
        end do
        value = value / scale
    end function relative_product_deviation

    pure subroutine build_angular_grids(n_theta, n_zeta, theta, zeta)
        integer, intent(in) :: n_theta, n_zeta
        real(dp), allocatable, intent(out) :: theta(:), zeta(:)
        integer :: i

        allocate (theta(n_theta), zeta(n_zeta))
        do i = 1, n_theta
            theta(i) = real(i - 1, dp) / real(n_theta, dp)
        end do
        do i = 1, n_zeta
            zeta(i) = real(i - 1, dp) / real(n_zeta, dp)
        end do
    end subroutine build_angular_grids

    pure subroutine allocate_result(ns, result)
        integer, intent(in) :: ns
        type(field_profile_identity_result_t), intent(out) :: result

        allocate (result%s(ns), result%toroidal_flux_deviation(ns))
        allocate (result%field_toroidal_flux_slope(ns))
        allocate (result%field_poloidal_flux_slope(ns))
        allocate (result%field_covariant_theta(ns))
        allocate (result%field_covariant_zeta(ns))
        allocate (result%poloidal_flux_deviation(ns))
        allocate (result%covariant_theta_deviation(ns))
        allocate (result%covariant_zeta_deviation(ns))
        allocate (result%iota_flux_deviation(ns))
        allocate (result%ampere_theta_deviation(ns))
        allocate (result%ampere_zeta_deviation(ns))
        allocate (result%exported_jacobian_deviation(ns))
        allocate (result%general_force_balance_deviation(ns))
    end subroutine allocate_result

end module field_profile_identities
