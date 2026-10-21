module perpendicular_kinetic_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ideal_mhd_energy_generated, only: generated_perpendicular_matrix, &
        generated_perpendicular_energy, generated_perpendicular_jvp, &
        generated_perpendicular_vjp
    use ideal_mhd_batch, only: mhd_mass_fp64, mhd_mass_fp32, mhd_mass_mixed
    implicit none
    private

    public :: benchmark_perpendicular_kinetic_energy
    public :: perpendicular_kinetic_energy
    public :: perpendicular_kinetic_matrix

    public :: perpendicular_kinetic_energy_jvp, perpendicular_kinetic_energy_vjp
    integer, parameter, public :: perpendicular_kinetic_input_count = 7

    public :: perpendicular_kinetic_energy_batch
    interface perpendicular_kinetic_energy_batch
        procedure mhd_mass_fp64, mhd_mass_fp32, mhd_mass_mixed
    end interface

contains

    pure subroutine perpendicular_kinetic_matrix(signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, density_kg_m3, mass)
        real(dp), intent(in) :: signed_sqrtg, bmag, grad_s2
        real(dp), intent(in) :: signed_sigma_tilde, density_kg_m3
        real(dp), intent(out) :: mass(2, 2)
        call generated_perpendicular_matrix(signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, density_kg_m3, mass)
    end subroutine perpendicular_kinetic_matrix

    pure function perpendicular_kinetic_energy(signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, density_kg_m3, displacement) result(energy)
        real(dp), intent(in) :: signed_sqrtg, bmag, grad_s2
        real(dp), intent(in) :: signed_sigma_tilde, density_kg_m3
        real(dp), intent(in) :: displacement(2)
        real(dp) :: energy, values(1)

        call generated_perpendicular_energy(signed_sqrtg, bmag, grad_s2, &
            signed_sigma_tilde, density_kg_m3, displacement(1), displacement(2), &
            values)
        energy = values(1)
    end function perpendicular_kinetic_energy

    pure function benchmark_perpendicular_kinetic_energy(active) result(energy) &
            bind(c, name="gliss_benchmark_perpendicular_kinetic_energy")
        real(c_double), intent(in), value :: active
        real(c_double) :: energy
        real(dp) :: displacement(2)

        displacement(1) = 0.4_dp + 0.02_dp * active
        displacement(2) = -0.3_dp + 0.01_dp * active
        energy = perpendicular_kinetic_energy(-1.1_dp + 0.01_dp * active, &
            1.4_dp + 0.04_dp * active, 1.3_dp + 0.02_dp * active, &
            0.2_dp - 0.01_dp * active, 2.3_dp + 0.05_dp * active, &
            displacement)
    end function benchmark_perpendicular_kinetic_energy

    ! Packed inputs follow the corresponding scalar energy/density argument order.
    ! Domain: nonzero signed Jacobian, positive B, grad_s2 and density.
    pure subroutine perpendicular_kinetic_energy_jvp(inputs, tangent, energy, derivative)
        real(dp), intent(in) :: inputs(7), tangent(7)
        real(dp), intent(out) :: energy, derivative
        real(dp) :: values(2)

        call generated_perpendicular_jvp(inputs(1), inputs(2), inputs(3), inputs(4), &
            inputs(5), inputs(6), inputs(7), tangent(1), tangent(2), tangent(3), &
            tangent(4), tangent(5), tangent(6), tangent(7), values)
        energy = values(1)
        derivative = values(2)
    end subroutine perpendicular_kinetic_energy_jvp

    pure subroutine perpendicular_kinetic_energy_vjp(inputs, cotangent, energy, gradient)
        real(dp), intent(in) :: inputs(7), cotangent
        real(dp), intent(out) :: energy, gradient(7)
        real(dp) :: values(8)

        call generated_perpendicular_vjp(inputs(1), inputs(2), inputs(3), inputs(4), &
            inputs(5), inputs(6), inputs(7), cotangent, values)
        energy = values(1)
        gradient = values(2:)
    end subroutine perpendicular_kinetic_energy_vjp

end module perpendicular_kinetic_kernel
