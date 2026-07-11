module perpendicular_kinetic_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: benchmark_perpendicular_kinetic_energy
    public :: perpendicular_kinetic_energy
    public :: perpendicular_kinetic_matrix

contains

    pure subroutine perpendicular_kinetic_matrix(signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, density_kg_m3, mass)
        real(dp), intent(in) :: signed_sqrtg, bmag, grad_s2
        real(dp), intent(in) :: signed_sigma_tilde, density_kg_m3
        real(dp), intent(out) :: mass(2, 2)
        real(dp) :: coefficients(2, 2)

        coefficients = 0.0_dp
        coefficients(1, 1) = 1.0_dp / sqrt(grad_s2)
        coefficients(2, 1) = signed_sigma_tilde / sqrt(grad_s2)
        coefficients(2, 2) = sqrt(grad_s2) / bmag
        mass = density_kg_m3 * abs(signed_sqrtg) &
            * matmul(transpose(coefficients), coefficients)
    end subroutine perpendicular_kinetic_matrix

    pure function perpendicular_kinetic_energy(signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, density_kg_m3, displacement) result(energy)
        real(dp), intent(in) :: signed_sqrtg, bmag, grad_s2
        real(dp), intent(in) :: signed_sigma_tilde, density_kg_m3
        real(dp), intent(in) :: displacement(2)
        real(dp) :: energy, mass(2, 2)

        call perpendicular_kinetic_matrix(signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, density_kg_m3, mass)
        energy = 0.5_dp * dot_product(displacement, &
            matmul(mass, displacement))
    end function perpendicular_kinetic_energy

    pure function benchmark_perpendicular_kinetic_energy(active) result(energy) &
            bind(c, name="gliss_benchmark_perpendicular_kinetic_energy")
        real(c_double), intent(in), value :: active
        real(c_double) :: energy
        real(dp) :: displacement(2)

        displacement = [0.4_dp + 0.02_dp * active, &
            -0.3_dp + 0.01_dp * active]
        energy = perpendicular_kinetic_energy(-1.1_dp + 0.01_dp * active, &
            1.4_dp + 0.04_dp * active, 1.3_dp + 0.02_dp * active, &
            0.2_dp - 0.01_dp * active, 2.3_dp + 0.05_dp * active, &
            displacement)
    end function benchmark_perpendicular_kinetic_energy

end module perpendicular_kinetic_kernel
