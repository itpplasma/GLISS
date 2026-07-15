module physical_mass_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: benchmark_physical_mass_energy
    public :: physical_mass_energy
    public :: physical_mass_matrix

contains

    pure subroutine physical_mass_matrix(flux_t_slope, flux_p_slope, &
            current_i, current_j, signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, &
            beta_tilde, density_kg_m3, mass)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: current_i, current_j, signed_sqrtg
        real(dp), intent(in) :: bmag, grad_s2, signed_sigma_tilde
        real(dp), intent(in) :: beta_tilde, density_kg_m3
        real(dp), intent(out) :: mass(3, 3)
        real(dp) :: coefficients(3, 3), flux_norm_squared, grad_s

        flux_norm_squared = flux_t_slope**2 + flux_p_slope**2
        grad_s = sqrt(grad_s2)
        coefficients = 0.0_dp
        coefficients(1, 1) = 1.0_dp / grad_s
        coefficients(2, 1) = signed_sigma_tilde / grad_s
        coefficients(2, 2) = grad_s / bmag
        coefficients(3, 1) = beta_tilde / bmag
        coefficients(3, 2) = -(current_i * flux_p_slope &
            - current_j * flux_t_slope) / (bmag * flux_norm_squared)
        coefficients(3, 3) = bmag / flux_norm_squared
        mass = density_kg_m3 * abs(signed_sqrtg) &
            * matmul(transpose(coefficients), coefficients)
    end subroutine physical_mass_matrix

    pure function physical_mass_energy(flux_t_slope, flux_p_slope, &
            current_i, current_j, signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, &
            beta_tilde, density_kg_m3, displacement) result(energy)
        real(dp), intent(in) :: flux_t_slope, flux_p_slope
        real(dp), intent(in) :: current_i, current_j, signed_sqrtg
        real(dp), intent(in) :: bmag, grad_s2, signed_sigma_tilde
        real(dp), intent(in) :: beta_tilde, density_kg_m3
        real(dp), intent(in) :: displacement(3)
        real(dp) :: energy, mass(3, 3), image(3)

        call physical_mass_matrix(flux_t_slope, flux_p_slope, current_i, &
            current_j, signed_sqrtg, bmag, grad_s2, signed_sigma_tilde, &
            beta_tilde, density_kg_m3, mass)
        image = matmul(mass, displacement)
        energy = 0.5_dp * dot_product(displacement, image)
    end function physical_mass_energy

    pure function benchmark_physical_mass_energy(active) result(energy) &
            bind(c, name="gvstab_benchmark_physical_mass_energy")
        real(c_double), intent(in), value :: active
        real(c_double) :: energy
        real(dp) :: displacement(3)

        displacement(1) = 0.4_dp + 0.02_dp * active
        displacement(2) = -0.3_dp + 0.01_dp * active
        displacement(3) = 0.2_dp - 0.03_dp * active
        energy = physical_mass_energy(1.2_dp + 0.05_dp * active, &
            0.7_dp - 0.03_dp * active, 0.8_dp + 0.02_dp * active, &
            0.6_dp - 0.01_dp * active, -1.1_dp + 0.01_dp * active, &
            1.4_dp + 0.04_dp * active, 1.3_dp + 0.02_dp * active, &
            0.2_dp - 0.01_dp * active, -0.15_dp + 0.03_dp * active, &
            2.3_dp + 0.05_dp * active, displacement)
    end function benchmark_physical_mass_energy

end module physical_mass_kernel
