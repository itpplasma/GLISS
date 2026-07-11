module two_component_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: bending_component_value
    public :: compression_component_value
    public :: shear_component_value
    public :: two_component_components
    public :: two_component_component_values
    public :: two_component_density

contains

    pure function magnetic_field_derivative(flux_t_slope, flux_p_slope, &
            sqrtg, angular_theta, angular_zeta) result(derivative)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope, sqrtg
        real(dp), intent(in) :: angular_theta, angular_zeta
        real(dp) :: derivative

        derivative = (flux_p_slope * angular_theta &
            + flux_t_slope * angular_zeta) / sqrtg
    end function magnetic_field_derivative

    pure function bending_component_value(flux_t_slope, flux_p_slope, &
            sqrtg, grad_s2, xi_s_theta, xi_s_zeta) result(component)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope, sqrtg, grad_s2
        real(dp), intent(in) :: xi_s_theta, xi_s_zeta
        real(dp) :: component

        component = magnetic_field_derivative(flux_t_slope, flux_p_slope, &
            sqrtg, xi_s_theta, xi_s_zeta) / sqrt(grad_s2)
    end function bending_component_value

    pure function shear_component_value(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, sqrtg, bmag, grad_s2, j_dot_b, &
            sigma_tilde, xi_s, xi_s_theta, xi_s_zeta, eta_theta, eta_zeta) &
            result(component)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve, sqrtg, bmag
        real(dp), intent(in) :: grad_s2, j_dot_b, sigma_tilde, xi_s
        real(dp), intent(in) :: xi_s_theta, xi_s_zeta, eta_theta, eta_zeta
        real(dp) :: component, bgrad_xi, bgrad_eta

        bgrad_xi = magnetic_field_derivative(flux_t_slope, flux_p_slope, &
            sqrtg, xi_s_theta, xi_s_zeta)
        bgrad_eta = magnetic_field_derivative(flux_t_slope, flux_p_slope, &
            sqrtg, eta_theta, eta_zeta)
        component = -(sqrt(grad_s2) / (bmag * sqrtg)) &
            * (sqrtg * bgrad_eta - (flux_t_slope * flux_p_curve &
            - flux_t_curve * flux_p_slope) * xi_s &
            + j_dot_b * sqrtg * xi_s / grad_s2 &
            + sigma_tilde * bmag * sqrtg * bgrad_xi / grad_s2)
    end function shear_component_value

    pure function compression_component_value(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, current_i, current_j, sqrtg, bmag, &
            pressure_slope, beta_tilde, xi_s, xi_s_s, xi_s_theta, &
            xi_s_zeta, eta_theta, eta_zeta) result(component)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve
        real(dp), intent(in) :: current_i, current_j, sqrtg, bmag
        real(dp), intent(in) :: pressure_slope, beta_tilde, xi_s, xi_s_s
        real(dp), intent(in) :: xi_s_theta, xi_s_zeta, eta_theta, eta_zeta
        real(dp) :: component, bgrad_xi

        bgrad_xi = magnetic_field_derivative(flux_t_slope, flux_p_slope, &
            sqrtg, xi_s_theta, xi_s_zeta)
        component = (current_j * eta_zeta - current_i * eta_theta &
            - (flux_t_slope * current_i + flux_p_slope * current_j) &
            * xi_s_s - (current_j * flux_p_curve &
            + current_i * flux_t_curve) * xi_s &
            - pressure_slope * sqrtg * xi_s &
            + beta_tilde * sqrtg * bgrad_xi) / (bmag * sqrtg)
    end function compression_component_value

    pure function two_component_component_values(flux_t_slope, &
            flux_p_slope, flux_t_curve, flux_p_curve, current_i, current_j, &
            sqrtg, bmag, grad_s2, j_dot_b, pressure_slope, sigma_tilde, &
            beta_tilde, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, &
            eta_zeta) result(components)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve
        real(dp), intent(in) :: current_i, current_j, sqrtg, bmag
        real(dp), intent(in) :: grad_s2, j_dot_b, pressure_slope
        real(dp), intent(in) :: sigma_tilde, beta_tilde
        real(dp), intent(in) :: xi_s, xi_s_s, xi_s_theta, xi_s_zeta
        real(dp), intent(in) :: eta_theta, eta_zeta
        real(dp) :: components(3)

        components(1) = bending_component_value(flux_t_slope, flux_p_slope, &
            sqrtg, grad_s2, xi_s_theta, xi_s_zeta)
        components(2) = shear_component_value(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, sqrtg, bmag, grad_s2, j_dot_b, &
            sigma_tilde, xi_s, xi_s_theta, xi_s_zeta, eta_theta, eta_zeta)
        components(3) = compression_component_value(flux_t_slope, &
            flux_p_slope, flux_t_curve, flux_p_curve, current_i, current_j, &
            sqrtg, bmag, pressure_slope, beta_tilde, xi_s, xi_s_s, &
            xi_s_theta, xi_s_zeta, eta_theta, eta_zeta)
    end function two_component_component_values

    pure subroutine two_component_components(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, current_i, current_j, sqrtg, &
            bmag, grad_s2, j_dot_b, pressure_slope, sigma_tilde, &
            beta_tilde, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, &
            eta_zeta, c_bending, c_shear, c_compression)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve
        real(dp), intent(in) :: current_i, current_j, sqrtg, bmag
        real(dp), intent(in) :: grad_s2, j_dot_b, pressure_slope
        real(dp), intent(in) :: sigma_tilde, beta_tilde
        real(dp), intent(in) :: xi_s, xi_s_s, xi_s_theta, xi_s_zeta
        real(dp), intent(in) :: eta_theta, eta_zeta
        real(dp), intent(out) :: c_bending, c_shear, c_compression
        real(dp) :: components(3)

        components = two_component_component_values(flux_t_slope, &
            flux_p_slope, flux_t_curve, flux_p_curve, current_i, current_j, &
            sqrtg, bmag, grad_s2, j_dot_b, pressure_slope, sigma_tilde, &
            beta_tilde, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, &
            eta_zeta)
        c_bending = components(1)
        c_shear = components(2)
        c_compression = components(3)
    end subroutine two_component_components

    pure subroutine two_component_density(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, current_i, current_j, sqrtg, &
            bmag, grad_s2, j_dot_b, pressure_slope, sigma_tilde, &
            beta_tilde, drive_a, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, &
            eta_theta, eta_zeta, density)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve
        real(dp), intent(in) :: current_i, current_j, sqrtg, bmag
        real(dp), intent(in) :: grad_s2, j_dot_b, pressure_slope
        real(dp), intent(in) :: sigma_tilde, beta_tilde, drive_a
        real(dp), intent(in) :: xi_s, xi_s_s, xi_s_theta, xi_s_zeta
        real(dp), intent(in) :: eta_theta, eta_zeta
        real(dp), intent(out) :: density
        real(dp) :: components(3)

        components = two_component_component_values(flux_t_slope, &
            flux_p_slope, flux_t_curve, flux_p_curve, current_i, current_j, &
            sqrtg, bmag, grad_s2, j_dot_b, pressure_slope, sigma_tilde, &
            beta_tilde, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, &
            eta_zeta)
        density = (components(1)**2 + components(2)**2 + components(3)**2 &
            - drive_a * xi_s**2) * abs(sqrtg)
    end subroutine two_component_density

end module two_component_kernel
