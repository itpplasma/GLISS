module three_component_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use physical_constants, only: vacuum_permeability
    use two_component_kernel, only: bending_component_value, &
        compression_component_value, shear_component_value
    implicit none
    private

    public :: benchmark_three_component_energy
    public :: compressible_divergence
    public :: compressible_divergence_value
    public :: three_component_density
    public :: three_component_density_value

contains

    pure function compressible_divergence_value(flux_t_slope, flux_p_slope, &
            signed_sqrtg, sqrtg_xi_s_radial, sqrtg_eta_theta, &
            sqrtg_eta_zeta, mu_theta, mu_zeta) result(divergence)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope, signed_sqrtg
        real(dp), intent(in) :: sqrtg_xi_s_radial, sqrtg_eta_theta
        real(dp), intent(in) :: sqrtg_eta_zeta, mu_theta, mu_zeta
        real(dp) :: divergence, flux_norm_squared

        flux_norm_squared = flux_t_slope**2 + flux_p_slope**2
        divergence = sqrtg_xi_s_radial / signed_sqrtg &
            + (flux_t_slope * sqrtg_eta_theta &
            - flux_p_slope * sqrtg_eta_zeta &
            + flux_p_slope * mu_theta + flux_t_slope * mu_zeta) &
            / (signed_sqrtg * flux_norm_squared)
    end function compressible_divergence_value

    pure subroutine compressible_divergence(flux_t_slope, flux_p_slope, &
            signed_sqrtg, sqrtg_xi_s_radial, sqrtg_eta_theta, &
            sqrtg_eta_zeta, mu_theta, mu_zeta, divergence)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope, signed_sqrtg
        real(dp), intent(in) :: sqrtg_xi_s_radial, sqrtg_eta_theta
        real(dp), intent(in) :: sqrtg_eta_zeta, mu_theta, mu_zeta
        real(dp), intent(out) :: divergence
        divergence = compressible_divergence_value(flux_t_slope, &
            flux_p_slope, signed_sqrtg, sqrtg_xi_s_radial, &
            sqrtg_eta_theta, sqrtg_eta_zeta, mu_theta, mu_zeta)
    end subroutine compressible_divergence

    pure function three_component_density_value(flux_t_slope, &
            flux_p_slope, flux_t_curve, flux_p_curve, current_i, current_j, &
            signed_sqrtg, bmag, grad_s2, j_dot_b, pressure_slope, &
            signed_sigma_tilde, beta_tilde, drive_a, gamma_pressure_pa, &
            xi_s, xi_s_radial, xi_s_theta, xi_s_zeta, eta_theta, eta_zeta, &
            sqrtg_xi_s_radial, sqrtg_eta_theta, sqrtg_eta_zeta, mu_theta, &
            mu_zeta) result(density)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve
        real(dp), intent(in) :: current_i, current_j, signed_sqrtg, bmag
        real(dp), intent(in) :: grad_s2, j_dot_b, pressure_slope
        real(dp), intent(in) :: signed_sigma_tilde, beta_tilde, drive_a
        real(dp), intent(in) :: gamma_pressure_pa, xi_s, xi_s_radial
        real(dp), intent(in) :: xi_s_theta, xi_s_zeta
        real(dp), intent(in) :: eta_theta, eta_zeta
        real(dp), intent(in) :: sqrtg_xi_s_radial, sqrtg_eta_theta
        real(dp), intent(in) :: sqrtg_eta_zeta, mu_theta, mu_zeta
        real(dp) :: density, c_bending, c_shear, c_compression, divergence

        c_bending = bending_component_value(flux_t_slope, flux_p_slope, &
            signed_sqrtg, grad_s2, xi_s_theta, xi_s_zeta)
        c_shear = shear_component_value(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, signed_sqrtg, bmag, grad_s2, &
            j_dot_b, signed_sigma_tilde, xi_s, xi_s_theta, xi_s_zeta, &
            eta_theta, eta_zeta)
        c_compression = compression_component_value(flux_t_slope, &
            flux_p_slope, flux_t_curve, flux_p_curve, current_i, current_j, &
            signed_sqrtg, bmag, pressure_slope, beta_tilde, xi_s, &
            xi_s_radial, xi_s_theta, xi_s_zeta, eta_theta, eta_zeta)
        divergence = compressible_divergence_value(flux_t_slope, &
            flux_p_slope, signed_sqrtg, sqrtg_xi_s_radial, &
            sqrtg_eta_theta, sqrtg_eta_zeta, mu_theta, mu_zeta)
        density = ((c_bending**2 + c_shear**2 + c_compression**2 &
            - drive_a * xi_s**2) &
            / vacuum_permeability + gamma_pressure_pa * divergence**2) &
            * abs(signed_sqrtg)
    end function three_component_density_value

    pure subroutine three_component_density(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, current_i, current_j, signed_sqrtg, &
            bmag, grad_s2, j_dot_b, pressure_slope, signed_sigma_tilde, &
            beta_tilde, drive_a, gamma_pressure_pa, xi_s, xi_s_radial, &
            xi_s_theta, xi_s_zeta, eta_theta, eta_zeta, &
            sqrtg_xi_s_radial, sqrtg_eta_theta, sqrtg_eta_zeta, &
            mu_theta, mu_zeta, density)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: flux_t_curve, flux_p_curve
        real(dp), intent(in) :: current_i, current_j, signed_sqrtg, bmag
        real(dp), intent(in) :: grad_s2, j_dot_b, pressure_slope
        real(dp), intent(in) :: signed_sigma_tilde, beta_tilde, drive_a
        real(dp), intent(in) :: gamma_pressure_pa, xi_s, xi_s_radial
        real(dp), intent(in) :: xi_s_theta, xi_s_zeta
        real(dp), intent(in) :: eta_theta, eta_zeta
        real(dp), intent(in) :: sqrtg_xi_s_radial, sqrtg_eta_theta
        real(dp), intent(in) :: sqrtg_eta_zeta, mu_theta, mu_zeta
        real(dp), intent(out) :: density
        density = three_component_density_value(flux_t_slope, flux_p_slope, &
            flux_t_curve, flux_p_curve, current_i, current_j, signed_sqrtg, &
            bmag, grad_s2, j_dot_b, pressure_slope, signed_sigma_tilde, &
            beta_tilde, drive_a, gamma_pressure_pa, xi_s, xi_s_radial, &
            xi_s_theta, xi_s_zeta, eta_theta, eta_zeta, &
            sqrtg_xi_s_radial, sqrtg_eta_theta, sqrtg_eta_zeta, mu_theta, &
            mu_zeta)
    end subroutine three_component_density

    pure function benchmark_three_component_energy(active) result(energy) &
            bind(c, name="gvstab_benchmark_three_component_energy")
        real(c_double), intent(in), value :: active
        real(c_double) :: energy

        energy = three_component_density_value(1.2_dp + 0.01_dp * active, &
            0.7_dp - 0.02_dp * active, 0.04_dp + 0.01_dp * active, &
            -0.03_dp + 0.01_dp * active, 0.8_dp + 0.02_dp * active, &
            0.6_dp - 0.01_dp * active, -1.1_dp + 0.01_dp * active, &
            1.4_dp + 0.02_dp * active, 1.3_dp + 0.01_dp * active, &
            0.2_dp - 0.01_dp * active, -0.15_dp + 0.01_dp * active, &
            0.1_dp - 0.01_dp * active, -0.12_dp + 0.02_dp * active, &
            0.05_dp + 0.01_dp * active, 0.9_dp + 0.03_dp * active, &
            0.3_dp + 0.01_dp * active, -0.2_dp + 0.02_dp * active, &
            0.1_dp - 0.01_dp * active, 0.04_dp + 0.02_dp * active, &
            -0.03_dp + 0.01_dp * active, 0.07_dp - 0.01_dp * active, &
            -0.05_dp + 0.01_dp * active, 0.08_dp + 0.01_dp * active, &
            -0.06_dp + 0.02_dp * active, 0.03_dp - 0.01_dp * active, &
            -0.02_dp + 0.01_dp * active)
    end function benchmark_three_component_energy

end module three_component_kernel
