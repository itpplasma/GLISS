module enzyme_physical_mass_family_fixture
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phase_assembly_policy, only: phase_assembly_transformed
    use physical_mass_family_assembly, only: &
        assemble_physical_family_mass_fixed_layout
    use radial_space_policy, only: radial_space_config_t
    implicit none
    private

    public :: physical_family_mass_energy

contains

    function physical_family_mass_energy(active) result(energy) &
            bind(c, name="physical_family_mass_energy")
        real(c_double), value :: active
        real(c_double) :: energy
        integer, parameter :: trial_m(2) = [1, 2]
        integer, parameter :: trial_n(2) = [1, 4]
        integer, parameter :: parity(2) = [1, 2]
        integer, parameter :: element_to_global(8, 3) = reshape([ &
            0, 0, 1, 2, 5, 6, 11, 12, &
            1, 2, 3, 4, 7, 8, 13, 14, &
            3, 4, 0, 0, 9, 10, 15, 16], [8, 3])
        real(dp), parameter :: stored_power(2) = [0.0_dp, 0.0_dp]
        type(radial_space_config_t) :: radial_space
        real(dp) :: fields(2, 2, 13, 3), density_kg_m3(3), mass(16, 16)
        real(dp) :: weight
        integer :: i, j, surface, row, column, info

        fields = 0.0_dp
        do surface = 1, 3
            do j = 1, 2
                do i = 1, 2
                    call set_active_fields(fields(i, j, :, surface), active, &
                        i, j, surface)
                end do
            end do
            density_kg_m3(surface) = 2.0_dp + 0.2_dp * surface &
                + active * 0.03_dp * surface
        end do
        call assemble_physical_family_mass_fixed_layout(fields, density_kg_m3, &
            trial_m, trial_n, parity, stored_power, 3, radial_space, &
            1.0_dp / 3.0_dp, phase_assembly_transformed, element_to_global, &
            mass, info)
        energy = 0.0_c_double
        if (info /= 0) return
        do column = 1, size(mass, 2)
            do row = 1, size(mass, 1)
                weight = real(row + 2 * column, dp) &
                    / real(3 * size(mass, 1), dp)
                energy = energy + weight * mass(row, column)
            end do
        end do
    end function physical_family_mass_energy

    pure subroutine set_active_fields(fields, active, i, j, surface)
        real(dp), intent(out) :: fields(:)
        real(c_double), intent(in) :: active
        integer, intent(in) :: i, j, surface

        fields = 0.0_dp
        fields(1) = 1.2_dp + 0.01_dp * surface &
            + active * 0.002_dp * (i + j)
        fields(2) = 0.7_dp + 0.01_dp * j - active * 0.001_dp * i
        fields(5) = 0.8_dp + active * 0.002_dp * i
        fields(6) = 0.6_dp - active * 0.001_dp * j
        fields(7) = -1.1_dp + active * 0.001_dp * (i + j)
        fields(8) = 1.4_dp + active * 0.002_dp * surface
        fields(9) = 1.3_dp + active * 0.001_dp * i
        fields(12) = 0.2_dp - active * 0.001_dp * j
        fields(13) = -0.15_dp + active * 0.002_dp * i
    end subroutine set_active_fields

end module enzyme_physical_mass_family_fixture
