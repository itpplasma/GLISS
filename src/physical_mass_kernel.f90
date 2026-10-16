module physical_mass_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gliss_physical_mass_generated, only: generated_physical_mass_matrix, &
        generated_physical_mass_energy, generated_physical_mass_jvp, &
        generated_physical_mass_vjp
    implicit none
    private

    public :: benchmark_physical_mass_energy
    public :: physical_mass_energy
    public :: physical_mass_matrix
    public :: physical_mass_energy_jvp, physical_mass_energy_vjp
    integer, parameter, public :: physical_mass_input_count = 13

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

        call generated_physical_mass_matrix(flux_t_slope, flux_p_slope, &
            current_i, current_j, signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, beta_tilde, density_kg_m3, mass)
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
        real(dp) :: energy

        call generated_physical_mass_energy(flux_t_slope, flux_p_slope, &
            current_i, current_j, signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, beta_tilde, density_kg_m3, displacement(1), &
            displacement(2), displacement(3), energy)
    end function physical_mass_energy

    ! Inputs and tangent order: FT', FP', I, J, signed sqrt(g), |B|,
    ! |grad(s)|^2, signed sigma_tilde, beta_tilde, rho, xi^s, eta, mu.
    ! Units and conventions match physical_mass_energy. Differentiation is
    ! local to J /= 0, B > 0, grad_s2 > 0, FT'^2+FP'^2 > 0, rho > 0.
    pure subroutine physical_mass_energy_jvp(inputs, tangent, energy, derivative)
        real(dp), intent(in) :: inputs(physical_mass_input_count)
        real(dp), intent(in) :: tangent(physical_mass_input_count)
        real(dp), intent(out) :: energy, derivative

        call generated_physical_mass_jvp(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), inputs(9), &
            inputs(10), inputs(11), inputs(12), inputs(13), tangent(1), &
            tangent(2), tangent(3), tangent(4), tangent(5), tangent(6), &
            tangent(7), tangent(8), tangent(9), tangent(10), tangent(11), &
            tangent(12), tangent(13), energy, derivative)
    end subroutine physical_mass_energy_jvp

    ! Same input order and open domain as physical_mass_energy_jvp.
    pure subroutine physical_mass_energy_vjp(inputs, cotangent, energy, gradient)
        real(dp), intent(in) :: inputs(physical_mass_input_count), cotangent
        real(dp), intent(out) :: energy, gradient(physical_mass_input_count)

        call generated_physical_mass_vjp(inputs(1), inputs(2), inputs(3), &
            inputs(4), inputs(5), inputs(6), inputs(7), inputs(8), inputs(9), &
            inputs(10), inputs(11), inputs(12), inputs(13), cotangent, &
            energy, gradient)
    end subroutine physical_mass_energy_vjp

    pure function benchmark_physical_mass_energy(active) result(energy) &
            bind(c, name="gvstab_benchmark_physical_mass_energy")
        real(c_double), intent(in), value :: active
        real(c_double) :: energy
        real(dp), parameter :: displacement_offset(3) = [0.4_dp, -0.3_dp, &
            0.2_dp]
        real(dp), parameter :: displacement_slope(3) = [0.02_dp, 0.01_dp, &
            -0.03_dp]
        real(dp) :: displacement(3)
        integer :: i

        do i = 1, 3
            displacement(i) = displacement_offset(i) &
                + displacement_slope(i) * active
        end do
        energy = physical_mass_energy(1.2_dp + 0.05_dp * active, &
            0.7_dp - 0.03_dp * active, 0.8_dp + 0.02_dp * active, &
            0.6_dp - 0.01_dp * active, -1.1_dp + 0.01_dp * active, &
            1.4_dp + 0.04_dp * active, 1.3_dp + 0.02_dp * active, &
            0.2_dp - 0.01_dp * active, -0.15_dp + 0.03_dp * active, &
            2.3_dp + 0.05_dp * active, displacement)
    end function benchmark_physical_mass_energy

end module physical_mass_kernel
