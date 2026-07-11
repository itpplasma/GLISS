module enzyme_physical_mass_assembly_fixture
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phase_assembly_policy, only: phase_assembly_transformed
    use physical_mass_assembly, only: &
        assemble_physical_mass_surface_resolved
    use radial_space_policy, only: radial_space_config_t
    implicit none
    private

    public :: transformed_physical_mass_energy

contains

    function transformed_physical_mass_energy(active) result(energy) &
            bind(c, name="transformed_physical_mass_energy")
        real(c_double), value :: active
        real(c_double) :: energy
        integer, parameter :: trial_m(3) = [1, 2, 1]
        integer, parameter :: trial_n(3) = [1, 4, -2]
        integer, parameter :: parity(3) = [1, 2, 1]
        real(dp), parameter :: stored_power(3) = [0.0_dp, 0.0_dp, 0.0_dp]
        type(radial_space_config_t) :: radial_space
        real(dp) :: fields(2, 2, 13), mass(12, 12), density_kg_m3
        real(dp) :: weight
        integer :: i, j, row, column, info

        fields = 0.0_dp
        do j = 1, 2
            do i = 1, 2
                fields(i, j, 1) = 1.2_dp + 0.01_dp * i &
                    + active * 0.002_dp * (i + j)
                fields(i, j, 2) = 0.7_dp + 0.01_dp * j &
                    - active * 0.001_dp * (2 * i + j)
                fields(i, j, 5) = 0.8_dp + active * 0.002_dp * i
                fields(i, j, 6) = 0.6_dp - active * 0.001_dp * j
                fields(i, j, 7) = -1.1_dp + active * 0.001_dp * (i + j)
                fields(i, j, 8) = 1.4_dp + active * 0.002_dp * j
                fields(i, j, 9) = 1.3_dp + active * 0.001_dp * i
                fields(i, j, 12) = 0.2_dp - active * 0.001_dp * j
                fields(i, j, 13) = -0.15_dp + active * 0.002_dp * i
            end do
        end do
        density_kg_m3 = 2.3_dp + active * 0.04_dp
        call assemble_physical_mass_surface_resolved(fields, density_kg_m3, &
            trial_m, trial_n, parity, stored_power, 3, radial_space, &
            0.375_dp, 0.25_dp, phase_assembly_transformed, mass, info)
        energy = 0.0_c_double
        if (info /= 0) return
        do column = 1, size(mass, 2)
            do row = 1, size(mass, 1)
                weight = real(row + 2 * column, dp) &
                    / real(3 * size(mass, 1), dp)
                energy = energy + weight * mass(row, column)
            end do
        end do
    end function transformed_physical_mass_energy

end module enzyme_physical_mass_assembly_fixture
