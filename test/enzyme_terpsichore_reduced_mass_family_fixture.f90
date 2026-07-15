module enzyme_terpsichore_reduced_mass_family_fixture
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use terpsichore_reduced_mass_family_assembly, only: &
        assemble_terpsichore_reduced_family_mass_fixed_layout
    implicit none
    private

    public :: terpsichore_reduced_family_mass_energy

contains

    function terpsichore_reduced_family_mass_energy(active) result(energy) &
            bind(c, name="terpsichore_reduced_family_mass_energy")
        real(c_double), value :: active
        real(c_double) :: energy
        integer, parameter :: element_to_global(6, 3) = reshape([ &
            0, 0, 1, 2, 5, 6, &
            1, 2, 3, 4, 7, 8, &
            3, 4, 0, 0, 9, 10], [6, 3])
        real(dp) :: signed_bjac(4, 3), flux_t_slope(3)
        real(dp) :: normal_phase(2, 4, 3), tangential_phase(2, 4, 3)
        real(dp) :: radial_factor(2, 3), radial_weight(3), mass(10, 10)
        real(dp) :: weight
        integer :: column, info, interval, point, row

        do interval = 1, 3
            flux_t_slope(interval) = 1.1_dp + 0.05_dp * interval &
                + 0.01_dp * active
            radial_weight(interval) = 0.4_dp + 0.1_dp * interval &
                - 0.02_dp * active
            radial_factor(1, interval) = 0.8_dp + 0.01_dp * interval &
                + 0.02_dp * active
            radial_factor(2, interval) = 1.1_dp - 0.01_dp * active
            do point = 1, 4
                signed_bjac(point, interval) = -(1.0_dp + 0.03_dp * point &
                    + 0.02_dp * interval + 0.01_dp * active * point)
                normal_phase(1, point, interval) = &
                    0.1_dp * point + 0.01_dp * active
                normal_phase(2, point, interval) = &
                    0.2_dp * interval - 0.02_dp * active
                tangential_phase(1, point, interval) = &
                    0.15_dp * interval + 0.02_dp * active
                tangential_phase(2, point, interval) = &
                    0.12_dp * point - 0.01_dp * active
            end do
        end do
        call assemble_terpsichore_reduced_family_mass_fixed_layout( &
            signed_bjac, flux_t_slope, normal_phase, tangential_phase, &
            radial_factor, radial_weight, element_to_global, mass, info)
        energy = 0.0_c_double
        if (info /= 0) return
        do column = 1, size(mass, 2)
            do row = 1, size(mass, 1)
                weight = real(row + 2 * column, dp) &
                    / real(3 * size(mass, 1), dp)
                energy = energy + weight * mass(row, column)
            end do
        end do
    end function terpsichore_reduced_family_mass_energy

end module enzyme_terpsichore_reduced_mass_family_fixture
