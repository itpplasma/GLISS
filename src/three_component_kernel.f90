module three_component_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use physical_constants, only: vacuum_permeability
    use ideal_mhd_energy_generated, only: generated_divergence, generated_three_energy, &
        generated_three_jvp, generated_three_vjp
    use ideal_mhd_batch, only: mhd_three_fp64, mhd_three_fp32, mhd_three_mixed
    implicit none
    private

    public :: benchmark_three_component_energy
    public :: compressible_divergence
    public :: compressible_divergence_value
    public :: three_component_density
    public :: three_component_density_value

    public :: three_component_energy_jvp, three_component_energy_vjp
    integer, parameter, public :: three_component_input_count = 26

    public :: three_component_energy_batch
    interface three_component_energy_batch
        procedure mhd_three_fp64, mhd_three_fp32, mhd_three_mixed
    end interface

contains

    pure function compressible_divergence_value(flux_t_slope, flux_p_slope, &
            signed_sqrtg, sqrtg_xi_s_radial, sqrtg_eta_theta, &
            sqrtg_eta_zeta, mu_theta, mu_zeta) result(divergence)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope, signed_sqrtg
        real(dp), intent(in) :: sqrtg_xi_s_radial, sqrtg_eta_theta
        real(dp), intent(in) :: sqrtg_eta_zeta, mu_theta, mu_zeta
        real(dp) :: divergence, values(1)

        call generated_divergence(flux_t_slope, flux_p_slope, signed_sqrtg, &
            sqrtg_xi_s_radial, sqrtg_eta_theta, sqrtg_eta_zeta, mu_theta, mu_zeta, &
            values)
        divergence = values(1)
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
        real(dp) :: density, values(1)

        call generated_three_energy(flux_t_slope, flux_p_slope, flux_t_curve, &
            flux_p_curve, current_i, current_j, signed_sqrtg, bmag, grad_s2, &
            j_dot_b, pressure_slope, signed_sigma_tilde, beta_tilde, drive_a, &
            gamma_pressure_pa, xi_s, xi_s_radial, xi_s_theta, xi_s_zeta, eta_theta, &
            eta_zeta, sqrtg_xi_s_radial, sqrtg_eta_theta, sqrtg_eta_zeta, mu_theta, &
            mu_zeta, vacuum_permeability, values)
        density = values(1)
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

    ! Packed inputs follow the corresponding scalar energy/density argument order.
    ! Domain: nonzero signed Jacobian, positive B and grad_s2; flux norm nonzero.
    pure subroutine three_component_energy_jvp(inputs, tangent, energy, derivative)
        real(dp), intent(in) :: inputs(26), tangent(26)
        real(dp), intent(out) :: energy, derivative
        real(dp) :: values(2)

        call generated_three_jvp(inputs(1), inputs(2), inputs(3), inputs(4), &
            inputs(5), inputs(6), inputs(7), inputs(8), inputs(9), inputs(10), &
            inputs(11), inputs(12), inputs(13), inputs(14), inputs(15), inputs(16), &
            inputs(17), inputs(18), inputs(19), inputs(20), inputs(21), inputs(22), &
            inputs(23), inputs(24), inputs(25), inputs(26), vacuum_permeability, &
            tangent(1), tangent(2), tangent(3), tangent(4), tangent(5), tangent(6), &
            tangent(7), tangent(8), tangent(9), tangent(10), tangent(11), &
            tangent(12), tangent(13), tangent(14), tangent(15), tangent(16), &
            tangent(17), tangent(18), tangent(19), tangent(20), tangent(21), &
            tangent(22), tangent(23), tangent(24), tangent(25), tangent(26), values)
        energy = values(1)
        derivative = values(2)
    end subroutine three_component_energy_jvp

    pure subroutine three_component_energy_vjp(inputs, cotangent, energy, gradient)
        real(dp), intent(in) :: inputs(26), cotangent
        real(dp), intent(out) :: energy, gradient(26)
        real(dp) :: values(27)

        call generated_three_vjp(inputs(1), inputs(2), inputs(3), inputs(4), &
            inputs(5), inputs(6), inputs(7), inputs(8), inputs(9), inputs(10), &
            inputs(11), inputs(12), inputs(13), inputs(14), inputs(15), inputs(16), &
            inputs(17), inputs(18), inputs(19), inputs(20), inputs(21), inputs(22), &
            inputs(23), inputs(24), inputs(25), inputs(26), vacuum_permeability, &
            cotangent, values)
        energy = values(1)
        gradient = values(2:)
    end subroutine three_component_energy_vjp

end module three_component_kernel
