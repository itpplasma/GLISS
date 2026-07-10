module two_component_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: two_component_components
    public :: two_component_density

contains

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
        real(dp) :: bgrad_xi, bgrad_eta

        bgrad_xi = (flux_p_slope * xi_s_theta &
            + flux_t_slope * xi_s_zeta) / sqrtg
        bgrad_eta = (flux_p_slope * eta_theta &
            + flux_t_slope * eta_zeta) / sqrtg
        c_bending = bgrad_xi / sqrt(grad_s2)
        c_shear = -(sqrt(grad_s2) / (bmag * sqrtg)) * (sqrtg * bgrad_eta &
            - (flux_t_slope * flux_p_curve &
            - flux_t_curve * flux_p_slope) * xi_s &
            + j_dot_b * sqrtg * xi_s / grad_s2 &
            + sigma_tilde * bmag * sqrtg * bgrad_xi / grad_s2)
        c_compression = (1.0_dp / (bmag * sqrtg)) * (current_j * eta_zeta &
            - current_i * eta_theta &
            - (flux_t_slope * current_i + flux_p_slope * current_j) &
            * xi_s_s &
            - (current_j * flux_p_curve + current_i * flux_t_curve) &
            * xi_s &
            - pressure_slope * sqrtg * xi_s &
            + beta_tilde * sqrtg * bgrad_xi)
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
        real(dp) :: c_bending, c_shear, c_compression

        call two_component_components(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, current_i, current_j, sqrtg, &
            bmag, grad_s2, j_dot_b, pressure_slope, sigma_tilde, &
            beta_tilde, xi_s, xi_s_s, xi_s_theta, xi_s_zeta, eta_theta, &
            eta_zeta, c_bending, c_shear, c_compression)
        density = (c_bending**2 + c_shear**2 + c_compression**2 &
            - drive_a * xi_s**2) * abs(sqrtg)
    end subroutine two_component_density

end module two_component_kernel
